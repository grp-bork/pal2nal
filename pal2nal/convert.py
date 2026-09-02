"""The peptide-to-codon matcher: pal2nal.v14.pl sub pn2codon, lines 1159-1919.

The approach is unchanged from v14. Each aligned residue contributes a
regex fragment matching the codons that could encode it, the fragments are
concatenated into one pattern, and that pattern is searched (unanchored,
case-insensitively) against the DNA. If the whole pattern matches, every
codon is correct by construction. If it does not, the peptide is cut into
anchors of ten residues; an anchor whose own pattern is found somewhere in
the DNA keeps its specific pattern, and one that is not found is relaxed to
wildcards. The mixed pattern is matched again and, this time, each codon is
checked individually so mismatches can be reported.

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

import re
from dataclasses import dataclass, field
from typing import Final

from .codontables import P2C

#: Residue characters v14 accepts as amino acids (pal2nal.pl:1762).
_AA_CHARS: Final = re.compile(r"[ACDEFGHIKLMNPQRSTVWY_*XU]")
_DIGIT: Final = re.compile(r"[0-9]")
_GAP: Final = re.compile(r"[-.]")

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


def _patterns(codontable: int) -> dict[str, str]:
    table = dict(P2C[codontable])
    table.setdefault("U", _SELENOCYSTEINE)
    return table


def _fragment(aa: str, p2c: dict[str, str]) -> str | None:
    """The regex fragment for one aligned character, or None for a gap."""
    if _AA_CHARS.fullmatch(aa):
        return p2c.get(aa, p2c["X"])
    if _DIGIT.fullmatch(aa):
        return "." * int(aa)
    if _GAP.fullmatch(aa):
        return None
    return p2c["X"]


def _build_pattern(pep: str, p2c: dict[str, str]) -> tuple[str, list[str]]:
    """pal2nal.pl:1757-1783. Returns the pattern and, separately, a report
    for every aligned character that is not a residue."""
    parts: list[str] = []
    messages: list[str] = []
    seen_letter = False
    for i, aa in enumerate(pep):
        pos = i + 1
        if _AA_CHARS.fullmatch(aa):
            if not seen_letter and aa == "M":
                # the initiating Met, which many tables spell more broadly
                parts.append(p2c["B"])
            else:
                parts.append(p2c.get(aa, p2c["X"]))
            seen_letter = True
        elif _DIGIT.fullmatch(aa):
            parts.append("." * int(aa))
        elif _GAP.fullmatch(aa):
            continue
        else:
            messages.append(f"pepAlnPos {pos}: {aa} unknown AA type. Taken as 'X'")
            parts.append(p2c["X"])
            seen_letter = True
    return "".join(parts), messages


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
        if not _GAP.fullmatch(aa):
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
    pattern, unknown = _build_pattern(pep, p2c)

    # exact match: every codon is right by construction, nothing to check
    if pattern:
        m = re.search(pattern, nuc, re.IGNORECASE)
        if m:
            return CodonResult(OK, m.group(0), unknown=unknown)

    # fallback: relax the anchors that cannot be found on their own
    whole: list[str] = []
    for index, anchor in enumerate(_anchors(pep)):
        specific: list[str] = []
        relaxed: list[str] = []
        first_anchor = index == 0
        seen_letter = False
        for aa in anchor:
            if _AA_CHARS.fullmatch(aa):
                if first_anchor and not seen_letter and aa == "M":
                    specific.append(p2c["B"])
                else:
                    specific.append(p2c.get(aa, p2c["X"]))
                relaxed.append(p2c["X"])
                seen_letter = True
            elif _DIGIT.fullmatch(aa):
                specific.append("." * int(aa))
                relaxed.append("." * int(aa))
            elif _GAP.fullmatch(aa):
                continue
            else:
                specific.append(p2c["X"])
                relaxed.append(p2c["X"])
                seen_letter = True
        joined = "".join(specific)
        whole.append(joined if joined and re.search(joined, nuc, re.IGNORECASE) else "".join(relaxed))

    combined = "".join(whole)
    if not combined:
        return CodonResult(NO_MATCH, unknown=unknown)
    m = re.search(combined, nuc, re.IGNORECASE)
    if not m:
        return CodonResult(NO_MATCH, unknown=unknown)

    codon = m.group(0)
    return CodonResult(MISMATCH, codon, _check_codons(pep, codon, p2c), unknown)


def _check_codons(pep: str, codon: str, p2c: dict[str, str]) -> list[str]:
    """pal2nal.pl:1878-1899 -- walk the match and verify each codon."""
    messages: list[str] = []
    pos_in_codon = 0
    residue_count = 0
    for i, aa in enumerate(pep):
        pos = i + 1
        if _DIGIT.fullmatch(aa):
            pos_in_codon += int(aa)
            continue
        if _GAP.fullmatch(aa):
            continue
        residue_count += 1
        tmpcodon = codon[pos_in_codon : pos_in_codon + 3]
        pos_in_codon += 3
        expected = p2c["B"] if residue_count == 1 and aa == "M" else p2c.get(aa)
        if expected is None or not re.search(expected, tmpcodon, re.IGNORECASE):
            messages.append(f"pepAlnPos {pos}: {aa} does not correspond to {tmpcodon}")
    return messages
