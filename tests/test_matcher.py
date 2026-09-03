"""Properties the codon matcher relies on, and a differential fuzz test.

`convert.pn2codon` reduces an aligned peptide to one regex per residue and
concatenates them. Two things make that work, and neither is obvious from
reading a single table:

* every fragment matches *exactly three* characters, so the whole pattern
  is fixed-width. That is what lets the matcher reason about codon
  boundaries by arithmetic, and what makes "the leftmost match" a
  well-defined offset rather than an artefact of backtracking;
* the fragments are a tiny language -- literals, ``.``, ``|`` and groups,
  over ``ACGTU`` plus IUPAC ambiguity codes -- with no quantifier,
  backreference or anchor anywhere.

Both are asserted here for all thirty-one tables, so a table added later
cannot quietly break the assumption. `test_fragments_accept_only_codons`
pins the alphabet the fragments can ever see, which is the alphabet
`validate.check_nucleotide_alphabet` lets through.

The fuzz test at the end compares `pn2codon` against a deliberately naive
reference matcher written straight from the definition: try every offset,
check every codon. It is far too slow for real input and exists only to
have a second opinion that shares no code with the implementation.
"""

from __future__ import annotations

import random
import re

import pytest

from pal2nal.codontables import P2C, SUPPORTED
from pal2nal.convert import MISMATCH, NO_MATCH, OK, _patterns, pn2codon
from pal2nal.validate import AMBIGUITY, NUCLEOTIDES

#: everything validate.check_nucleotide_alphabet admits, in both cases
DNA_ALPHABET = NUCLEOTIDES + AMBIGUITY
ALL_CODONS = [
    a + b + c for a in DNA_ALPHABET for b in DNA_ALPHABET for c in DNA_ALPHABET
]


@pytest.mark.parametrize("tid", SUPPORTED)
def test_every_fragment_is_exactly_three_wide(tid: int) -> None:
    """The invariant the whole matcher rests on."""
    for aa, fragment in P2C[tid].items():
        lo, hi = re._parser.parse(fragment).getwidth()
        assert (lo, hi) == (3, 3), f"table {tid}, {aa!r}: {fragment!r} is not 3 wide"


@pytest.mark.parametrize("tid", SUPPORTED)
def test_fragments_use_no_regex_feature_beyond_alternation(tid: int) -> None:
    """No quantifiers, anchors, backreferences or character classes: the
    fragments must stay a language a non-regex matcher could also read."""
    allowed = set("().|ACGTURYSWKMBDHVN")
    for aa, fragment in P2C[tid].items():
        stray = set(fragment) - allowed
        assert not stray, f"table {tid}, {aa!r}: unexpected {sorted(stray)}"


@pytest.mark.parametrize("tid", SUPPORTED)
def test_every_fragment_matches_at_least_one_real_codon(tid: int) -> None:
    """A fragment that matches nothing would silently desynchronise every
    residue after it, which is exactly what v14's empty U pattern did."""
    for aa, fragment in P2C[tid].items():
        assert any(
            re.fullmatch(fragment, c, re.IGNORECASE) for c in ALL_CODONS
        ), f"table {tid}, {aa!r}: matches no codon at all"


def test_the_degenerate_position_really_is_degenerate() -> None:
    """``.`` is a fully degenerate position and X is a fully degenerate
    codon, so they match characters the DNA validator would have refused.
    Nothing outside the alphabet can reach the matcher, so this is safe --
    but any replacement matcher has to reproduce it rather than assume the
    fragments only ever see ACGTU and ambiguity codes."""
    assert re.fullmatch(P2C[1]["X"], "@@@")
    assert re.fullmatch(P2C[1]["A"], "GC@")


@pytest.mark.parametrize("tid", SUPPORTED)
def test_fragments_are_upper_case_only(tid: int) -> None:
    """Lower-case DNA is accepted and matched under re.IGNORECASE, which
    only behaves symmetrically because no fragment spells a base in lower
    case. A stray lower-case letter would still match under IGNORECASE
    here but would break any matcher that upper-cases the DNA instead."""
    for aa, fragment in P2C[tid].items():
        assert fragment == fragment.upper(), f"table {tid}, {aa!r}: {fragment!r}"


def test_every_fragment_matches_the_same_codons_in_either_case() -> None:
    """The symmetry itself, checked exhaustively for the universal code;
    test_fragments_are_upper_case_only carries it to the other thirty."""
    for aa, fragment in P2C[1].items():
        for codon in ALL_CODONS:
            assert bool(re.fullmatch(fragment, codon, re.IGNORECASE)) == bool(
                re.fullmatch(fragment, codon.lower(), re.IGNORECASE)
            ), f"{aa!r}: {codon} and {codon.lower()} differ"


# --------------------------------------------------------------------------
# differential fuzz


def reference_match(pep: str, nuc: str, tid: int) -> str | None:
    """The definition, written the slow obvious way: the leftmost offset at
    which every residue's fragment matches its codon in turn.

    Shares nothing with convert.py but the tables themselves.
    """
    p2c = P2C[tid]
    seleno = "((U|T)GA)"
    steps: list[tuple[int, str | None]] = []
    seen_letter = False
    for aa in pep:
        if aa in "-.":
            continue
        if aa.isdigit():
            steps.append((int(aa), None))
            continue
        if aa == "U" and "U" not in p2c:
            steps.append((3, seleno))
            seen_letter = True
            continue
        if aa in p2c or aa in "ACDEFGHIKLMNPQRSTVWY_*XU":
            frag = p2c["B"] if aa == "M" and not seen_letter else p2c.get(aa, p2c["X"])
        else:
            frag = p2c["X"]
        steps.append((3, frag))
        seen_letter = True

    width = sum(w for w, _ in steps)
    for start in range(len(nuc) - width + 1):
        pos = start
        for w, frag in steps:
            if frag is not None and not re.fullmatch(
                frag, nuc[pos : pos + w], re.IGNORECASE
            ):
                break
            pos += w
        else:
            return nuc[start : start + width]
    return None


def random_case(rng: random.Random) -> tuple[str, str, int]:
    """A peptide alignment row and some DNA, sometimes corresponding."""
    tid = rng.choice(SUPPORTED)
    p2c = P2C[tid]
    residues = [a for a in p2c if a not in "BX"] + list("-.")
    pep = "".join(rng.choice(residues) for _ in range(rng.randint(1, 25)))
    if all(a in "-." for a in pep):        # no residues at all: see below
        pep = "M" + pep
    if rng.random() < 0.15:
        pep = pep[:-1] + rng.choice("123456789")
    if rng.random() < 0.6:
        # DNA that really does encode the peptide, with a random flank
        parts = []
        for aa in pep:
            if aa in "-.":
                continue
            if aa.isdigit():
                parts.append("".join(rng.choice("ACGT") for _ in range(int(aa))))
                continue
            frag = p2c["X"] if aa == "X" else p2c.get(aa, p2c["X"])
            hits = [c for c in ALL_CODONS if re.fullmatch(frag, c, re.IGNORECASE)]
            parts.append(rng.choice(hits))
        flank = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 9)))
        nuc = flank + "".join(parts) + "".join(rng.choice("ACGT") for _ in range(6))
    else:
        nuc = "".join(rng.choice(DNA_ALPHABET) for _ in range(rng.randint(3, 90)))
    if rng.random() < 0.2:
        nuc = nuc.lower()
    return pep, nuc, tid


@pytest.mark.parametrize("seed", range(40))
def test_pn2codon_agrees_with_the_naive_matcher(seed: int) -> None:
    """On an exact match, pn2codon must return the codons the definition
    picks out. Where no exact match exists it drops into the anchor
    fallback, which is a reporting path rather than a matching one, so only
    the OK verdict is compared."""
    rng = random.Random(seed)
    for _ in range(15):
        pep, nuc, tid = random_case(rng)
        got = pn2codon(pep, nuc, tid)
        want = reference_match(pep, nuc, tid)
        if want is not None:
            assert got.result == OK, f"{pep!r} vs {nuc!r} (table {tid})"
            assert got.codonseq == want, f"{pep!r} vs {nuc!r} (table {tid})"
        else:
            assert got.result in (MISMATCH, NO_MATCH), f"{pep!r} vs {nuc!r} (table {tid})"


def test_the_fuzz_generator_reaches_both_verdicts() -> None:
    """Guard against the fuzz test quietly testing only one branch."""
    rng = random.Random(0)
    seen = set()
    for _ in range(150):
        pep, nuc, tid = random_case(rng)
        seen.add(pn2codon(pep, nuc, tid).result)
    assert OK in seen and (MISMATCH in seen or NO_MATCH in seen), seen


@pytest.mark.parametrize("pep", ["-", ".", "---", "-.-.-"])
def test_a_row_of_nothing_but_gaps_finds_no_match(pep: str) -> None:
    """An all-gap alignment row builds an empty pattern. It is excluded
    from the fuzz comparison because "the leftmost offset matching nothing"
    is offset 0 by the definition and NO_MATCH here; this pins which of the
    two the port does, so an optimisation cannot quietly swap them."""
    assert pn2codon(pep, "ATGGCT", 1).result == NO_MATCH


def test_selenocysteine_is_matched_and_v14_left_it_undefined() -> None:
    """U has no pattern in v14's tables; convert.py supplies TGA."""
    assert all("U" not in P2C[t] for t in SUPPORTED)
    assert _patterns(1)["U"]
    assert pn2codon("MU", "ATGTGA", 1).result == OK


# --------------------------------------------------------------------------
# the bit-parallel matcher against the regexes it replaces


@pytest.mark.parametrize("tid", SUPPORTED)
def test_parsed_fragments_accept_exactly_what_the_regex_accepts(tid: int) -> None:
    """codonmatch reads the tables itself instead of handing them to re.

    The tables are the specification and stay written as regexes, so the
    reading has to agree with a regex engine on every codon -- including
    the ones outside the DNA alphabet that a degenerate position lets
    through, and both cases of every letter.
    """
    from pal2nal.codonmatch import codon_matches

    probes = [c for c in ALL_CODONS]
    probes += [c.lower() for c in ALL_CODONS[::37]]
    probes += ["@@@", "GC@", "A", "", "ACGT", "gc.", "..."]
    for aa, fragment in P2C[tid].items():
        for codon in probes:
            assert codon_matches(fragment, codon) == bool(
                re.fullmatch(fragment, codon, re.IGNORECASE)
            ), f"table {tid}, {aa!r}: {codon!r}"


@pytest.mark.parametrize("seed", range(30))
def test_offsets_search_agrees_with_an_unanchored_regex(seed: int) -> None:
    """`Offsets.search` must return what `re.search` on the concatenated
    pattern returned: the same leftmost offset, or None."""
    from pal2nal.codonmatch import Offsets

    rng = random.Random(1000 + seed)
    for _ in range(20):
        tid = rng.choice(SUPPORTED)
        p2c = P2C[tid]
        keys = [a for a in p2c if a != "B"]
        steps: list[tuple[int, str | None]] = []
        for _ in range(rng.randint(1, 12)):
            if rng.random() < 0.12:
                width = rng.randint(1, 9)
                steps.append((width, None))
            else:
                steps.append((3, p2c[rng.choice(keys)]))
        nuc = "".join(rng.choice("ACGTURYN") for _ in range(rng.randint(0, 60)))
        if rng.random() < 0.3:
            nuc = nuc.lower()
        pattern = "".join(
            "." * width if frag is None else frag for width, frag in steps
        )
        want = re.search(pattern, nuc, re.IGNORECASE)
        got = Offsets(nuc).search(steps)
        assert got == (want.start() if want else None), f"{steps!r} in {nuc!r}"


def test_offsets_handles_an_empty_sequence() -> None:
    from pal2nal.codonmatch import Offsets

    empty = Offsets("")
    assert empty.length == 0
    assert empty.search([(3, P2C[1]["A"])]) is None


def test_a_fragment_the_parser_cannot_read_is_refused() -> None:
    """The parser trusts the fragment grammar; anything else must fail
    loudly rather than silently match the wrong codons."""
    from pal2nal.codonmatch import branches

    with pytest.raises(ValueError):
        branches("(GC.)(GC.)")      # six wide
    with pytest.raises(ValueError):
        branches("(GC.")            # unbalanced
