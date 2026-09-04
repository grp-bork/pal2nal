"""The -nogap filter must use the selected genetic code's stop codons.

v14 wrote the universal code's TAA/TAG/TGA into the filter literally and
ignored -codontable, so on every other code the filter was wrong in both
directions at once: it deleted columns that spell an amino acid, and it kept
the in-frame stops that -nogap exists to remove. The golden corpus pins that
behaviour for two probes; this checks the fix against every supported table,
codon by codon, using the same published references test_codontables.py
diffs the tables themselves against.
"""

from __future__ import annotations

import pytest
from test_codontables import CODONS, DATA, parse_ncbi, parse_wikipedia

from pal2nal.codontables import SUPPORTED
from pal2nal.output import Alignment, remove_gaps


@pytest.fixture(scope="module")
def published() -> dict[int, list[str]]:
    """{table id: its stop codons}, from NCBI's gc.prt and Wikipedia."""
    ncbi = parse_ncbi((DATA / "ncbi_gc.prt").read_text())
    wiki = parse_wikipedia((DATA / "wikipedia_gc_34_37.tsv").read_text())
    stops = {tid: list(t["stops"]) for tid, t in ncbi.items()}
    for tid, aa in wiki.items():
        stops[tid] = [c for c, a in aa.items() if a == "*"]
    return stops


def _filtered(codons: list[str], codontable: int) -> list[str]:
    """The codons -nogap keeps, given one sequence holding all of them."""
    aln = Alignment(ids=["s"])
    aln.codonaln = ["".join(codons)]
    aln.coloraln = list(aln.codonaln)
    aln.maskseq = "-" * len(aln.codonaln[0])
    kept = remove_gaps(aln, codontable).codonaln[0]
    return [kept[i : i + 3] for i in range(0, len(kept), 3)]


@pytest.mark.parametrize("tid", SUPPORTED)
def test_nogap_drops_exactly_the_tables_stop_codons(tid, published):
    assert tid in published, f"no published reference for table {tid}"
    expected = [c for c in CODONS if c not in published[tid]]
    assert _filtered(CODONS, tid) == expected


def test_nogap_defaults_to_the_universal_code():
    """cli.py passes the table through; the default must stay table 1."""
    assert _filtered(CODONS, 1) == [c for c in CODONS if c not in ("TAA", "TAG", "TGA")]
    aln = Alignment(ids=["s"])
    aln.codonaln = ["".join(CODONS)]
    aln.coloraln = list(aln.codonaln)
    aln.maskseq = "-" * len(aln.codonaln[0])
    assert remove_gaps(aln).codonaln == remove_gaps(aln, 1).codonaln


@pytest.mark.parametrize(
    "tid,codon,keep",
    [
        # the two directions the v14 filter got wrong, spelled out
        (6, "TAA", True),  # Gln under the ciliate code, not a stop
        (6, "TAG", True),  # Gln
        (6, "TGA", False),  # table 6's only stop
        (2, "TGA", True),  # Trp under the vertebrate mitochondrial code
        (2, "AGA", False),  # a stop v14's universal regex never matched
        (2, "AGG", False),
        (14, "TAA", True),  # Tyr under the alternative flatworm code
        (14, "TAG", False),
    ],
)
def test_known_disagreements_with_the_universal_code(tid, codon, keep):
    assert (codon in _filtered([codon], tid)) is keep


def test_a_gapped_column_still_goes(published):
    """The gap half of the filter is unchanged by the codon-table fix."""
    assert _filtered(["AT-", "ATG"], 1) == ["ATG"]
