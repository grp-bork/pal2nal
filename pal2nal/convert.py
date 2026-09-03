"""The peptide-to-codon matcher: pal2nal.v14.pl sub pn2codon, lines 1159-1919.

The approach is unchanged from v14. Each aligned residue contributes a
pattern matching the codons that could encode it, the patterns are
concatenated, and the result is searched (unanchored, case-insensitively)
against the DNA. If the whole thing matches, every codon is correct by
construction. If it does not, the peptide is cut into anchors of ten
residues; an anchor whose own pattern is found somewhere in the DNA keeps
its specific pattern, and one that is not found is relaxed to wildcards.
The mixed pattern is matched again and, this time, each codon is checked
individually so mismatches can be reported.

What is not v14's is who does the searching. The patterns are still the
transcribed genetic codes of `codontables.P2C`, written as regexes, but
`codonmatch` reads them directly instead of handing a peptide-length
pattern to `re` for every sequence. See that module for why; the answers
are the same ones `re.search` gave, which `tests/test_matcher.py` checks
fragment by fragment and offset by offset.

Divergences from v14, all agreed and listed in PORTING.md:

* a peptide short enough to produce a single anchor no longer loses that
  anchor to an empty loop range (v14 left the pattern undefined, so any
  alignment of ten residues or fewer that needed this path failed outright)
* U (selenocysteine) has a codon pattern instead of an empty one
* the fallback path uses the selected table's own initiation codons rather
  than a hardcoded set that matches none of the seventeen tables
* "." counts as a gap everywhere, as the documented contract always said
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from typing import Final

from .codonmatch import Offsets, codon_matches
from .codontables import P2C

#: Residue characters v14 accepts as amino acids (pal2nal.pl:1762).
_AA_CHARS: Final = frozenset("ACDEFGHIKLMNPQRSTVWY_*XU")
_DIGIT: Final = frozenset("0123456789")
_GAP: Final = frozenset("-.")

#: A width and the codon pattern that must match there, or None for a
#: frame-shift numeral, which consumes nucleotides without constraining
#: them. Concatenating these is what v14 did by joining pattern strings.
Step = tuple[int, "str | None"]

#: v14 lists U among the accepted residues but defines no pattern for it,
#: so it contributed nothing and desynchronised everything after it. TGA is
#: the selenocysteine codon.
_SELENOCYSTEINE: Final = "((U|T)GA)"

OK: Final = 1
MISMATCH: Final = 2
NO_MATCH: Final = -1


@dataclass
class CodonResult:
    result: int
    codonseq: str = ""
    #: codons that do not encode the aligned residue
    messages: list[str] = field(default_factory=list)
    #: aligned characters that are not residues at all, taken as X. Kept
    #: apart from the mismatches because they are an input problem, and
    #: v15 reports them whatever the output filters say
    unknown: list[str] = field(default_factory=list)


@cache
def _patterns(codontable: int) -> dict[str, str]:
    table = dict(P2C[codontable])
    table.setdefault("U", _SELENOCYSTEINE)
    return table


def _build_steps(pep: str, p2c: dict[str, str]) -> tuple[list[Step], list[str]]:
    """pal2nal.pl:1757-1783. Returns the steps and, separately, a report
    for every aligned character that is not a residue."""
    steps: list[Step] = []
    messages: list[str] = []
    seen_letter = False
    for i, aa in enumerate(pep):
        pos = i + 1
        if aa in _AA_CHARS:
            if not seen_letter and aa == "M":
                # the initiating Met, which many tables spell more broadly
                steps.append((3, p2c["B"]))
            else:
                steps.append((3, p2c.get(aa, p2c["X"])))
            seen_letter = True
        elif aa in _DIGIT:
            steps.append((int(aa), None))
        elif aa in _GAP:
            continue
        else:
            messages.append(f"pepAlnPos {pos}: {aa} unknown AA type. Taken as 'X'")
            steps.append((3, p2c["X"]))
            seen_letter = True
    return steps, messages


def _width(steps: list[Step]) -> int:
    """How many nucleotides the steps consume. Fixed, never a range: that
    is the property the whole matcher rests on."""
    return sum(width for width, _ in steps)


def _anchors(pep: str) -> list[str]:
    """pal2nal.pl:1807-1832 -- ten non-gap residues per anchor.

    v14 merged a short trailing anchor into its predecessor with a loop over
    `0..$#preanchor - 1`, which is empty when there is only one anchor; the
    anchor list came out empty and the match was made against an undefined
    pattern. A single anchor is now kept as-is.
    """
    pre: list[str] = []
    current: list[str] = []
    count = 0
    for i, aa in enumerate(pep):
        current.append(aa)
        if aa not in _GAP:
            count += 1
        if count == 10 or i == len(pep) - 1:
            pre.append("".join(current))
            current = []
            count = 0
    if not pre:
        return []
    if len(pre) > 1 and len(pre[-1]) < 10:
        merged = pre[:-2] + [pre[-2] + pre[-1]]
        return merged
    return pre


def pn2codon(pep: str, nuc: str, codontable: int) -> CodonResult:
    p2c = _patterns(codontable)
    steps, unknown = _build_steps(pep, p2c)
    # one pass over the DNA, shared by the exact match and every anchor
    offsets = Offsets(nuc)

    # exact match: every codon is right by construction, nothing to check
    if steps:
        start = offsets.search(steps)
        if start is not None:
            # sliced from the original, so the DNA's own case survives
            return CodonResult(OK, nuc[start : start + _width(steps)], unknown=unknown)

    # fallback: relax the anchors that cannot be found on their own
    whole: list[Step] = []
    for index, anchor in enumerate(_anchors(pep)):
        specific: list[Step] = []
        relaxed: list[Step] = []
        first_anchor = index == 0
        seen_letter = False
        for aa in anchor:
            if aa in _AA_CHARS:
                if first_anchor and not seen_letter and aa == "M":
                    specific.append((3, p2c["B"]))
                else:
                    specific.append((3, p2c.get(aa, p2c["X"])))
                relaxed.append((3, p2c["X"]))
                seen_letter = True
            elif aa in _DIGIT:
                specific.append((int(aa), None))
                relaxed.append((int(aa), None))
            elif aa in _GAP:
                continue
            else:
                specific.append((3, p2c["X"]))
                relaxed.append((3, p2c["X"]))
                seen_letter = True
        found = bool(specific) and offsets.search(specific) is not None
        whole.extend(specific if found else relaxed)

    if not whole:
        return CodonResult(NO_MATCH, unknown=unknown)
    start = offsets.search(whole)
    if start is None:
        return CodonResult(NO_MATCH, unknown=unknown)

    codon = nuc[start : start + _width(whole)]
    return CodonResult(MISMATCH, codon, _check_codons(pep, codon, p2c), unknown)


def _check_codons(pep: str, codon: str, p2c: dict[str, str]) -> list[str]:
    """pal2nal.pl:1878-1899 -- walk the match and verify each codon."""
    messages: list[str] = []
    pos_in_codon = 0
    residue_count = 0
    for i, aa in enumerate(pep):
        pos = i + 1
        if aa in _DIGIT:
            pos_in_codon += int(aa)
            continue
        if aa in _GAP:
            continue
        residue_count += 1
        tmpcodon = codon[pos_in_codon : pos_in_codon + 3]
        pos_in_codon += 3
        expected = p2c["B"] if residue_count == 1 and aa == "M" else p2c.get(aa)
        if expected is None or not codon_matches(expected, tmpcodon):
            messages.append(f"pepAlnPos {pos}: {aa} does not correspond to {tmpcodon}")
    return messages
