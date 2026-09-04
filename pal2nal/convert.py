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
#: v16, -partial only: some of the peptide was placed against the DNA and
#: the rest was gapped. Never returned unless the caller asks for it.
PARTIAL: Final = 3


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
    #: alignment columns -partial could not place, and so gapped. Columns
    #: rather than prose: a hole can span hundreds of residues, and one
    #: warning each would drown the report the way v14's per-codon
    #: mismatches do. cli.py renders them as a count and a range
    unplaced: list[int] = field(default_factory=list)
    #: runs of contiguous DNA the placement used. More than one means the
    #: DNA carried something the peptide does not -- an intron, typically
    segments: int = 0


@cache
def _patterns(codontable: int) -> dict[str, str]:
    table = dict(P2C[codontable])
    table.setdefault("U", _SELENOCYSTEINE)
    return table


def _has_letter(pep: str) -> bool:
    """Whether `pep` holds anything `_build_steps` would count as a letter.

    Frame-shift numerals and gaps do not count; an unknown character does,
    because it is taken as X. This is what lets `_chain` tell each anchor
    whether the initiating Met has already gone by.
    """
    return any(aa not in _GAP and aa not in _DIGIT for aa in pep)


def _build_steps(
    pep: str, p2c: dict[str, str], *, seen_letter: bool = False
) -> tuple[list[Step], list[str]]:
    """pal2nal.pl:1757-1783. Returns the steps and, separately, a report
    for every aligned character that is not a residue.

    `seen_letter` says whether a residue has already been seen *earlier in
    the peptide* than the slice passed in. It matters only for the
    initiating Met: `_chain` builds one anchor at a time, and without this
    every anchor beginning with M would be matched against the table's
    initiation codons instead of the Met codon. Under table 1 that would
    silently accept CTG and TTG as Met; under table 11, ATT, ATC, ATA and
    GTG as well.
    """
    steps: list[Step] = []
    messages: list[str] = []
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


def _chain(
    pep: str, nuc: str, offsets: Offsets, p2c: dict[str, str]
) -> tuple[str, list[int], int]:
    """Place each anchor on its own, left to right, and gap the rest.

    v14's fallback concatenates every anchor into one *fixed-width* pattern,
    so a single length discrepancy anywhere -- an intron, an indel, a CDS
    truncated at either end -- makes the whole search fail and the run
    abort. Placing the anchors independently, each at the leftmost offset
    that does not overlap the one before it, lets the DNA between two
    anchors be any length at all, which is precisely what those four
    pathologies need.

    Returns the codon string and the peptide indices it could not place.
    The codon string always spends exactly the width each aligned column
    expects -- three nucleotides per residue, `int(aa)` for a frame-shift
    numeral -- so `output.build` reinserts the alignment's gaps around it
    unchanged. That width is `_width(steps)`, never `3 * residues`: for
    "EK4QKNDTY" the true width is 28 and `3 * 9` is 27, and a nucleotide
    lost here would shift every codon downstream of it.
    """
    anchors = _anchors(pep)
    # (peptide offset, steps, width) per anchor, in peptide order
    plan: list[tuple[int, list[Step], int]] = []
    pep_at = 0
    seen_letter = False
    for anchor in anchors:
        steps, _ = _build_steps(anchor, p2c, seen_letter=seen_letter)
        plan.append((pep_at, steps, _width(steps)))
        pep_at += len(anchor)
        seen_letter = seen_letter or _has_letter(anchor)

    # pass 1: greedy left to right. On contiguous DNA the leftmost offset
    # at or after `pos` is `pos` itself, so a clean pair chains to exactly
    # the answer the unanchored search would have given.
    placed: list[int | None] = [None] * len(plan)
    pos = 0
    for i, (_, steps, width) in enumerate(plan):
        if not steps:
            continue
        start = offsets.search(steps, pos)
        if start is not None:
            placed[i] = start
            pos = start + width

    # pass 2: fill the runs pass 1 left behind. An anchor fails to place
    # for two very different reasons -- the DNA genuinely does not encode
    # it, or it does but with a substitution the pattern will not accept --
    # and the span its neighbours leave tells them apart.
    _fill_runs(plan, placed)

    pieces: list[str] = []
    gapped: list[int] = []
    segments = 0
    previous_end: int | None = None
    for i, (pep_at, steps, width) in enumerate(plan):
        at = placed[i]
        if at is None:
            pieces.append("-" * width)
            gapped.extend(
                pep_at + k
                for k, aa in enumerate(anchors[i])
                if aa not in _GAP and aa not in _DIGIT
            )
            previous_end = None
        else:
            pieces.append(nuc[at : at + width])
            if at != previous_end:
                segments += 1
            previous_end = at + width
    return "".join(pieces), gapped, segments


def _fill_runs(
    plan: list[tuple[int, list[Step], int]], placed: list[int | None]
) -> None:
    """Recover unplaced anchors whose position is pinned anyway.

    A run of unplaced anchors *between* two placed ones has its position
    determined from both sides: if the DNA they leave is exactly as wide as
    the run needs, those are its codons, whatever they encode. Filling
    there is what keeps -partial from being worse than the fallback it
    supplements -- an anchor holding one point substitution will not place,
    and without this it would be gapped rather than reported as a mismatch.

    Only interior runs. A run at either end is pinned on one side only, so
    the span proves nothing and a fill there would read 5' UTR as coding
    sequence -- which is exactly where the two truncation pathologies would
    trigger it. Verifying such a fill instead of trusting the span would
    not help: pass 1 searches from the same offset with the same patterns,
    so anything that verified there would already have been placed.
    """
    i = 0
    while i < len(plan):
        if placed[i] is not None:
            i += 1
            continue
        j = i
        while j < len(plan) and placed[j] is None:
            j += 1
        if i == 0 or j == len(plan):
            i = j                       # an end run: nothing pins it
            continue
        prev = placed[i - 1]
        assert prev is not None
        start = prev + plan[i - 1][2]
        need = sum(plan[k][2] for k in range(i, j))
        if placed[j] == start + need:   # an indel would not fit
            for k in range(i, j):
                placed[k] = start
                start += plan[k][2]
        i = j


def pn2codon(
    pep: str, nuc: str, codontable: int, *, partial: bool = False
) -> CodonResult:
    """Match `pep` against `nuc` and return its codons.

    `partial` is v16's -partial: where v14 would give up entirely, place
    each anchor on its own instead and gap what will not place. It is
    consulted only after both of v14's paths have failed, so every input
    that converts today converts to the same bytes either way.
    """
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

    start = offsets.search(whole) if whole else None
    if start is None:
        if not partial:
            return CodonResult(NO_MATCH, unknown=unknown)
        codon, gapped, segments = _chain(pep, nuc, offsets, p2c)
        messages = _check_codons(pep, codon, p2c, frozenset(gapped))
        return CodonResult(PARTIAL, codon, messages, unknown, gapped, segments)

    codon = nuc[start : start + _width(whole)]
    return CodonResult(MISMATCH, codon, _check_codons(pep, codon, p2c), unknown)


def _check_codons(
    pep: str, codon: str, p2c: dict[str, str], skip: frozenset[int] = frozenset()
) -> list[str]:
    """pal2nal.pl:1878-1899 -- walk the match and verify each codon.

    `skip` holds the columns -partial gapped, so a residue that has no
    codon at all is never reported as one whose codon is wrong. The cursor
    still advances across them: that arithmetic is the register the whole
    output depends on.
    """
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
        if i in skip:
            continue
        expected = p2c["B"] if residue_count == 1 and aa == "M" else p2c.get(aa)
        if expected is None or not codon_matches(expected, tmpcodon):
            messages.append(f"pepAlnPos {pos}: {aa} does not correspond to {tmpcodon}")
    return messages
