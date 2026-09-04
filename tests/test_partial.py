"""-partial: recovering a codon alignment from an imperfect pep/CDS pair.

The golden corpus pins what -partial *prints*. These are the properties it
cannot express: that the flag is inert on everything that already worked,
that the codon string keeps the alignment's register, and that the two
hazards in the placement code stay fixed.
"""

from __future__ import annotations

import io
import shlex
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from pal2nal.codonmatch import Offsets, codon_matches
from pal2nal.convert import (
    MISMATCH,
    NO_MATCH,
    OK,
    PARTIAL,
    _build_steps,
    _chain,
    _patterns,
    _width,
    pn2codon,
)
from pal2nal.output import build

TESTS = Path(__file__).resolve().parent
DATA = TESTS / "data"

#: one unambiguous codon per residue, so a construction below reads as the
#: protein it encodes
CODE = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT",
    "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT",
    "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
    "W": "TGG", "Y": "TAT",
}

#: every pathology that reaches v14's "inconsistency" error
PATHOLOGIES = ["intron", "indel", "trunc5", "trunc3", "partial_frameshift"]


def fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif line.strip():
            records[name] += line.strip()
    return records


def pair(stem: str) -> tuple[str, str]:
    """The single peptide and DNA of a one-sequence fixture."""
    pep = fasta(DATA / f"{stem}_pep.fasta")
    dna = fasta(DATA / f"{stem}.nuc")
    (sid,) = pep
    return pep[sid], dna[sid]


def residues(pep: str) -> int:
    return sum(1 for aa in pep if aa not in "-." and not aa.isdigit())


def run_cli(args: list[str]) -> tuple[str, str, int]:
    from pal2nal.cli import main

    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        status = main(args)
    return out.getvalue(), err.getvalue(), status


# --------------------------------------------------------------- the gate


@pytest.mark.parametrize("stem", PATHOLOGIES)
def test_partial_engages_only_where_the_run_would_abort(stem: str) -> None:
    """Every pathology aborts without the flag and converts with it."""
    pep, dna = pair(stem)
    assert pn2codon(pep, dna, 1).result == NO_MATCH
    assert pn2codon(pep, dna, 1, partial=True).result == PARTIAL


def test_an_unrelated_cds_is_left_to_the_old_fallback() -> None:
    """-partial does not rescue a DNA that never aborted in the first place.

    An unrelated CDS reaches v14's relaxed-wildcard fallback, which matches
    all-wildcards at offset 0 and "succeeds". That is a poor answer, but the
    evidence behind it is identical to the evidence behind a correct one --
    a short peptide carrying a single substitution matches the same way and
    is recovered properly. So the flag stays out of it, and the coverage
    line is what tells the user the result is worthless.
    """
    pep, dna = pair("unrelated")
    assert pn2codon(pep, dna, 1).result == MISMATCH
    assert pn2codon(pep, dna, 1, partial=True).result == MISMATCH
    _, err, status = run_cli(
        [str(DATA / "unrelated_pep.fasta"), str(DATA / "unrelated.nuc"), "-partial"]
    )
    assert status == 0
    assert "MISMATCH: seq1 1/60 codons verified" in err


@pytest.mark.parametrize("stem", ["simple", "mismatch"])
def test_the_ok_and_mismatch_paths_ignore_partial(stem: str) -> None:
    pep = fasta(DATA / f"{stem}_pep.fasta")
    dna = fasta(DATA / f"{stem}.nuc")
    for sid, seq in pep.items():
        plain = pn2codon(seq, dna[sid], 1)
        with_flag = pn2codon(seq, dna[sid], 1, partial=True)
        assert plain.result in (OK, MISMATCH)
        assert plain == with_flag


# -------------------------------------------------------------- recovery


def test_an_intron_recovers_the_exact_cds() -> None:
    """The codons either side of an intron are the true CDS, byte for byte."""
    pep, introned = pair("intron")
    # trunc3 is cds[:-45] and trunc5 is cds[45:], so together they are the CDS
    truth = pair("trunc3")[1][:45] + pair("trunc5")[1]
    res = pn2codon(pep, introned, 1, partial=True)
    assert res.codonseq == truth
    assert res.unplaced == []
    assert res.segments == 2, "the placement had to step over the intron"
    assert "N" not in res.codonseq.upper(), "no intron nucleotide leaked through"


def test_a_truncated_cds_gaps_the_end_it_does_not_cover() -> None:
    """Rather than fabricating codons out of whatever flanks the match."""
    for stem, first, last in (("trunc5", 0, 20), ("trunc3", 40, 60)):
        pep, dna = pair(stem)
        res = pn2codon(pep, dna, 1, partial=True)
        assert res.unplaced == list(range(first, last)), stem
        assert res.messages == [], f"{stem}: a hole is not a mismatch"


def test_unplaced_residues_are_holes_not_mismatches() -> None:
    """The two are reported apart, and together they account for everything."""
    for stem in PATHOLOGIES:
        pep, dna = pair(stem)
        res = pn2codon(pep, dna, 1, partial=True)
        placed = residues(pep) - len(res.unplaced)
        assert 0 <= placed <= residues(pep), stem
        # a column is a hole or a mismatch, never both
        holes = set(res.unplaced)
        for message in res.messages:
            col = int(message.split()[1].rstrip(":")) - 1
            assert col not in holes, f"{stem}: column {col} reported twice"


# ------------------------------------------------- the register invariant


@pytest.mark.parametrize("stem", PATHOLOGIES)
def test_the_codon_string_always_spends_the_full_width(stem: str) -> None:
    """What output.build's cursor arithmetic assumes, for every pathology."""
    pep, dna = pair(stem)
    expected = _width(_build_steps(pep, _patterns(1))[0])
    assert len(pn2codon(pep, dna, 1, partial=True).codonseq) == expected


def test_an_unplaced_frameshift_numeral_keeps_its_own_width() -> None:
    """The gap filler is `"-" * _width`, never three dashes per residue.

    A frame-shift numeral consumes `int(aa)` nucleotides. Here the peptide's
    first anchor holds a `4` and cannot be placed, so it is gapped; the
    second anchor can be, and follows it. Filling three dashes per residue
    spends 30 where the anchor needs 31, and every codon after the hole
    slides one nucleotide left -- the drift defect, one level up.
    """
    first, second = "M4KQLRTYWC", "AKQLRTYWCM"
    pep = first + second
    # the DNA encodes only the second anchor, so the first cannot be placed
    dna = "".join(CODE[aa] for aa in second)

    expected_width = 3 + 4 + 8 * 3          # M, the numeral, then 8 residues
    assert _width(_build_steps(first, _patterns(1))[0]) == expected_width == 31

    res = pn2codon(pep, dna, 1, partial=True)
    assert res.result == PARTIAL
    assert res.codonseq[:31] == "-" * 31, "the numeral's own width, not 30"
    assert res.codonseq[31:] == dna, "the placed anchor kept its register"

    # and the same thing seen through the output, where it would show up as
    # codons sitting under the wrong residues
    row = build([">x"], [pep], [res.codonseq], set(), "").codonaln[0]
    assert row.endswith(dna)


# ------------------------------------------------------- the initiating Met


def test_the_initiation_codon_applies_only_to_the_first_residue() -> None:
    """`_chain` builds one anchor at a time, so `seen_letter` is threaded.

    Under table 1 the initiation pattern `B` is `((U|T|C|Y|A)(U|T)G)`: it
    accepts CTG and TTG as well as ATG. `_build_steps` resets its own
    `seen_letter` per call, so calling it per anchor without threading would
    apply `B` to the first M of *every* anchor and let a mid-sequence Met
    match either of two Leu codons. Under table 11 it is five codons wide.
    """
    p2c = _patterns(1)
    assert codon_matches(p2c["B"], "CTG"), "B is the looser pattern"
    assert not codon_matches(p2c["M"], "CTG"), "M is not"

    opening, _ = _build_steps("MAKQ", p2c, seen_letter=False)
    later, _ = _build_steps("MAKQ", p2c, seen_letter=True)
    assert opening[0][1] == p2c["B"]
    assert later[0][1] == p2c["M"]

    # and end to end: an M at the head of the second anchor, encoded CTG.
    # Threaded, that anchor cannot be placed; unthreaded, B accepts CTG and
    # ten wrong codons are emitted as though they were verified.
    first, second = "MAKQLRTYWC", "MKQLRTYWCA"
    pep = first + second
    dna = (
        "TTG" + "".join(CODE[aa] for aa in first[1:])      # only B accepts TTG
        + "N" * 20                                         # an intron
        + "CTG" + "".join(CODE[aa] for aa in second[1:])
    )
    codonseq, gapped, _segments = _chain(pep, dna, Offsets(dna), p2c)
    assert gapped == list(range(10, 20)), (
        "the second anchor's M was matched with the initiation codons"
    )
    assert codonseq[:30] == "TTG" + "".join(CODE[aa] for aa in first[1:])
    assert codonseq[30:] == "-" * 30


# ------------------------------------------------- the point of the flag


def test_one_bad_pair_no_longer_stops_the_others() -> None:
    """Two good sequences convert identically whether or not a third,
    unmatchable one is in the alignment with them."""
    both = [str(DATA / "onebad_pep.fasta"), str(DATA / "onebad.nuc")]
    good = [
        str(DATA / "onebad_good_only_pep.fasta"),
        str(DATA / "onebad_good_only.nuc"),
    ]

    # without the flag the bad pair takes the whole run down with it
    _, _, status = run_cli(both)
    assert status == 1

    with_bad, err, status = run_cli(both + ["-partial", "-output", "fasta"])
    assert status == 0
    assert "UNMATCHED: bad 0/60 residues placed" in err
    without_bad, _, _ = run_cli(good + ["-output", "fasta"])

    def records(text: str) -> dict[str, str]:
        out: dict[str, str] = {}
        name = ""
        for line in text.splitlines():
            if line.startswith(">"):
                name = line[1:]
                out[name] = ""
            elif line.strip():
                out[name] += line.strip()
        return out

    kept, alone = records(with_bad), records(without_bad)
    assert set(kept) == {"good1", "bad", "good2"}
    assert set(alone) == {"good1", "good2"}
    for sid in ("good1", "good2"):
        assert kept[sid] == alone[sid], f"{sid} changed because of the bad pair"
    assert set(kept["bad"]) == {"-"}, "the unmatchable row is all gaps"


# --------------------------------------------------- inert where unneeded


def corpus_cases() -> list[tuple[str, str]]:
    rows = []
    for line in (TESTS / "cases.tsv").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, _, args = line.partition("\t")
        if not name.startswith("partial_"):
            rows.append((name, args))
    return rows


RULE = "#" + "-" * 72 + "#"
_COVERAGE = ("#  MISMATCH: ", "#  PARTIAL: ", "#  UNMATCHED: ")


def split_report(text: str) -> tuple[list[str], list[str]]:
    """Separate the "#---#" message block from the alignment around it.

    Under -html the block goes to stdout inside the <pre> rather than to
    stderr, so a caller that wants "did the alignment change" has to take
    the block out of whichever stream it landed in. A CLUSTAL mask line is
    also all hashes, hence the two-space prefix rather than one hash.
    """
    body, report = [], []
    after_rule = False
    for line in text.splitlines():
        # the block writer closes with the rule and a blank line; when the
        # block appears only because -partial had something to say, that
        # blank line is part of it and not of the alignment
        if line.startswith("#  ") or line == RULE or (after_rule and not line):
            report.append(line)
        else:
            body.append(line)
        after_rule = line == RULE or (after_rule and not line)
    return body, report


@pytest.mark.parametrize("name,args", corpus_cases())
def test_partial_never_changes_an_alignment_that_already_converted(
    name: str, args: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The gate, checked across the corpus rather than argued for.

    -partial is consulted only after both of v14's matching paths have
    failed, so on every input that converts today the codon alignment and
    the exit status must be untouched. The message block may gain coverage
    lines -- that is the flag reporting how much of each sequence its DNA
    accounted for, and under -nomismatch, which suppresses the warnings, it
    can be the only thing in the block. Nothing already in the block may
    change.

    Skipped: cases that abort without the flag, which are exactly the ones
    -partial exists to change, and cases that print the usage text, which
    is selected by the argument *count*, so appending a token changes which
    branch runs before any of this is reached.
    """
    monkeypatch.chdir(TESTS)
    argv = shlex.split(args)
    base_out, base_err, base_status = run_cli(argv)
    if "Usage:  pal2nal.pl" in base_out + base_err:
        pytest.skip("usage text is selected by the argument count")
    if base_status != 0:
        pytest.skip("aborts without -partial: changing that is the feature")

    out, err, status = run_cli([*argv, "-partial"])
    assert status == base_status

    base_body, base_report = split_report(base_out + base_err)
    body, report = split_report(out + err)
    assert body == base_body, "the codon alignment changed"

    # under -nomismatch the base run has no block at all, so a coverage line
    # brings the rule and the "Input files:" header along with it
    scaffolding = (RULE, "", "#  Input files:  ")
    added = [line for line in report if line not in base_report]
    assert all(
        line.startswith(_COVERAGE) or line.startswith(scaffolding) for line in added
    ), f"the report gained something other than a coverage line: {added}"
    messages = [
        line
        for line in report
        if not line.startswith(_COVERAGE) and not line.startswith(scaffolding)
    ]
    assert messages == [
        line
        for line in base_report
        if not line.startswith(_COVERAGE) and not line.startswith(scaffolding)
    ], "an existing message changed"
