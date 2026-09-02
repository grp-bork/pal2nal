"""Check every genetic code in codontables.py against its published source.

The tables are regexes that match codons, so a defect in one is invisible
until someone converts an alignment that happens to use the affected codon.
This decodes each table back into a plain codon -> amino acid mapping and
compares it against the reference data vendored in ``tests/data/``.

This is the check that found the table 10 defect: pal2nal.v14.pl left TGA in
the stop pattern and out of Cys, so ``-codontable 10`` silently returned
universal-code results. Run against the Perl's own hashes it reports exactly
one wrong codon in seventeen tables, which is how the bug was pinned down;
``test_perl_v14_has_exactly_one_table_defect`` keeps that finding on record.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pal2nal.codontables import (
    ASSUMED_STARTS,
    MODIFIED_TABLES,
    NAMES,
    P2C,
    PROVENANCE,
    SUPPORTED,
)

DATA = Path(__file__).parent / "data"
LEGACY = Path(__file__).parent / "reference" / "pal2nal.v14.pl"

BASES = "TCAG"
#: NCBI's ncbieaa order: first base slowest, third base fastest.
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]

#: Tables v14 shipped, i.e. those whose patterns must not drift.
V14_TABLES = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)


def parse_ncbi(text: str) -> dict[int, dict[str, object]]:
    """Read NCBI's gc.prt into {id: {name, aa, starts, stops}}.

    A codon counts as a stop if it is ``*`` in ncbieaa (an unconditional
    stop) or in sncbieaa (a context-dependent one, as in tables 27/28/31).
    """
    tables: dict[int, dict[str, object]] = {}
    for block in re.findall(r"\{([^{}]*?)\}", text, re.DOTALL):
        ids = re.search(r"\bid\s+(\d+)", block)
        aas = re.search(r'ncbieaa\s+"([^"]+)"', block)
        starts = re.search(r'sncbieaa\s+"([^"]+)"', block)
        name = re.search(r'name\s+"([^"]+)"', block)
        if not (ids and aas and starts):
            continue
        aa, st = aas.group(1), starts.group(1)
        assert len(aa) == len(st) == 64, f"table {ids.group(1)} is not 64 codons"
        tables[int(ids.group(1))] = {
            "name": name.group(1) if name else "",
            "aa": dict(zip(CODONS, aa)),
            "starts": [c for c, s in zip(CODONS, st) if s == "M"],
            "stops": [c for c, (a, s) in zip(CODONS, zip(aa, st)) if a == "*" or s == "*"],
        }
    return tables


def parse_wikipedia(text: str) -> dict[int, dict[str, str]]:
    """Read the vendored codes 34-37 into {id: {codon: amino acid}}."""
    out = {}
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        tid, aa = line.split("\t")
        assert len(aa) == 64, f"table {tid} is not 64 codons"
        out[int(tid)] = dict(zip(CODONS, aa))
    return out


def decode(table: dict[str, str]) -> dict[str, set[str]]:
    """Map each codon to the residue keys whose pattern accepts it.

    ``B`` (initiation codons), ``X`` (anything) and ``_`` (a synonym of
    ``*``) overlap the real residues by design and are left out.
    """
    compiled = {aa: re.compile(p + r"\Z") for aa, p in table.items() if aa not in ("B", "X", "_")}
    return {c: {aa for aa, rx in compiled.items() if rx.match(c)} for c in CODONS}


def starts_of(table: dict[str, str]) -> list[str]:
    rx = re.compile(table["B"] + r"\Z")
    return [c for c in CODONS if rx.match(c)]


@pytest.fixture(scope="module")
def ncbi() -> dict[int, dict[str, object]]:
    return parse_ncbi((DATA / "ncbi_gc.prt").read_text())


@pytest.fixture(scope="module")
def wikipedia() -> dict[int, dict[str, str]]:
    return parse_wikipedia((DATA / "wikipedia_gc_34_37.tsv").read_text())


def reference(tid, ncbi, wikipedia) -> dict[str, str]:
    return ncbi[tid]["aa"] if tid in ncbi else wikipedia[tid]


# --------------------------------------------------------------------------
# the tables themselves
# --------------------------------------------------------------------------


def test_metadata_is_consistent():
    assert set(SUPPORTED) == set(P2C) == set(NAMES) == set(PROVENANCE)
    assert len(SUPPORTED) == len(set(SUPPORTED)), "duplicate id in SUPPORTED"
    assert list(SUPPORTED) == sorted(SUPPORTED), "SUPPORTED is not in ascending order"


@pytest.mark.parametrize("tid", SUPPORTED)
def test_every_codon_is_assigned_exactly_once(tid):
    """No codon may be unmatched, and none may match two residues.

    The exception is a context-dependent stop, which legitimately matches
    both its amino acid and ``*``.
    """
    ncbi_data = parse_ncbi((DATA / "ncbi_gc.prt").read_text())
    dual = set()
    if tid in ncbi_data:
        dual = {c for c in ncbi_data[tid]["stops"] if ncbi_data[tid]["aa"][c] != "*"}
    for codon, hits in decode(P2C[tid]).items():
        if codon in dual:
            assert hits == {ncbi_data[tid]["aa"][codon], "*"}, (
                f"table {tid} codon {codon} should match its residue and '*', got {sorted(hits)}"
            )
        else:
            assert len(hits) == 1, (
                f"table {tid} codon {codon} matches {sorted(hits)}, expected exactly one"
            )


@pytest.mark.parametrize("tid", SUPPORTED)
def test_amino_acids_match_the_published_table(tid, ncbi, wikipedia):
    want = reference(tid, ncbi, wikipedia)
    got = decode(P2C[tid])
    wrong = [
        (c, want[c], sorted(got[c])) for c in CODONS if want[c] not in got[c]
    ]
    assert not wrong, f"table {tid} ({NAMES[tid]}) disagrees with its source: {wrong}"


@pytest.mark.parametrize("tid", SUPPORTED)
def test_initiation_codons_match_the_published_table(tid, ncbi):
    """Every table's ``B`` pattern accepts exactly the published start codons.

    Tables 34-37 have none published, so they take the bacterial set; see
    ASSUMED_STARTS.
    """
    if tid in ASSUMED_STARTS:
        assert tid not in ncbi, f"table {tid} now has NCBI data; drop the assumed starts"
        assert starts_of(P2C[tid]) == ncbi[11]["starts"], (
            f"table {tid} should assume the bacterial start set until NCBI publishes one"
        )
    else:
        assert starts_of(P2C[tid]) == ncbi[tid]["starts"], (
            f"table {tid} initiation codons disagree with NCBI"
        )


def test_tables_34_to_37_are_not_ncbi_assignments(ncbi):
    """Guard the numbering caveat: if NCBI ever assigns these, revisit them."""
    assert max(ncbi) == 33, (
        "the vendored gc.prt now defines ids past 33; codes 34-37 use Wikipedia's "
        "numbering and may collide, so check PROVENANCE before updating the file"
    )


# --------------------------------------------------------------------------
# the v14 tables must not drift, and the one v14 defect stays on record
# --------------------------------------------------------------------------


#: Initiation codons v15 accepts that pal2nal.v14.pl does not. v14's lists
#: predate later NCBI revisions; tables 3 and 13 gained theirs in gc.prt 4.4
#: and 4.0, and table 11's omission has no NCBI changelog entry at all.
V14_START_CODON_ADDITIONS = {
    3: ["GTG"],
    11: ["ATA", "ATC", "ATT"],
    13: ["ATA"],
}


@pytest.mark.parametrize("tid", V14_TABLES)
def test_v14_tables_are_unchanged_except_where_recorded(tid, ncbi):
    """The ported tables still agree with v14 except where MODIFIED_TABLES says.

    Both halves of a table count: the residue a codon codes for, and whether
    it can initiate translation.
    """
    perl = perl_tables()
    ours, theirs = decode(P2C[tid]), decode(perl[tid])
    residues = [c for c in CODONS if ours[c] != theirs[c]]
    starts = sorted(set(starts_of(P2C[tid])) ^ set(starts_of(perl[tid])))
    if tid in MODIFIED_TABLES:
        assert residues or starts, f"table {tid} is listed in MODIFIED_TABLES but matches v14"
    else:
        assert not residues, (
            f"table {tid} drifted from pal2nal.v14.pl at {residues}; "
            "record it in MODIFIED_TABLES if the change is intended"
        )
        assert not starts, (
            f"table {tid} initiation codons drifted from pal2nal.v14.pl at {starts}; "
            "record it in MODIFIED_TABLES if the change is intended"
        )


@pytest.mark.parametrize("tid", V14_TABLES)
def test_start_codons_widened_exactly_where_recorded(tid):
    """v15 adds initiation codons to three tables and takes none away."""
    perl = perl_tables()
    ours, theirs = set(starts_of(P2C[tid])), set(starts_of(perl[tid]))
    assert not theirs - ours, f"table {tid} dropped v14 initiation codons {sorted(theirs - ours)}"
    assert sorted(ours - theirs) == V14_START_CODON_ADDITIONS.get(tid, []), (
        f"table {tid} initiation codons changed unexpectedly: added {sorted(ours - theirs)}"
    )


def perl_tables() -> dict[int, dict[str, str]]:
    """Pull the %p2c hashes straight out of pal2nal.v14.pl."""
    text = LEGACY.read_text()
    blocks: dict[int, dict[str, str]] = {}
    for m in re.finditer(r"\$codontable == (\d+)\)\s*\{(.*?)\n        \);", text, re.DOTALL):
        entries = dict(re.findall(r'"([^"]+)"\s*=>\s*"([^"]+)"', m.group(2)))
        if entries:
            blocks.setdefault(int(m.group(1)), entries)
    first = re.search(r"%p2c = \((.*?)\n        \);", text, re.DOTALL)
    blocks.setdefault(1, dict(re.findall(r'"([^"]+)"\s*=>\s*"([^"]+)"', first.group(1))))
    return blocks


def test_perl_v14_has_exactly_one_table_defect(ncbi):
    """How the table 10 bug was found, kept runnable.

    Decoding v14's own hashes against NCBI turns up one wrong codon in the
    seventeen tables it shipped: TGA in the Euplotid nuclear code, which v14
    leaves as a stop where NCBI has it coding for Cys.
    """
    perl = perl_tables()
    assert sorted(perl) == list(V14_TABLES)
    defects = []
    for tid, table in perl.items():
        got = decode(table)
        for codon in CODONS:
            want = ncbi[tid]["aa"][codon]
            if want not in got[codon]:
                defects.append((tid, codon, want, sorted(got[codon])))
    assert defects == [(10, "TGA", "C", ["*"])], defects
    # the residue defect is table 10's alone; the other entries in
    # MODIFIED_TABLES are the initiation-codon corrections, which are
    # omissions rather than wrong assignments.
    assert 10 in MODIFIED_TABLES
    assert set(MODIFIED_TABLES) - {10} == set(V14_START_CODON_ADDITIONS)
