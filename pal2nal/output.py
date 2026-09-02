"""Building and printing the codon alignment: pal2nal.v14.pl lines 665-1156.

`pn2codon` returns each sequence's codons with no gaps in them. This module
puts the gaps back by walking the peptide alignment column by column, then
applies the -blockonly/-nomismatch/-nogap filters and renders one of the
four output formats.

Column widths are set by the widest entry in each column, so a frame-shift
that consumes four nucleotides in one sequence widens that column for every
sequence and the rows stay in register.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Final, TextIO

_DIGIT: Final = re.compile(r"[0-9]")
_GAP: Final = re.compile(r"[-.]")
#: in-frame stop codons removed by -nogap (pal2nal.pl:857)
_STOP: Final = re.compile(r"(((U|T)A(A|G|R))|((T|U)GA))", re.IGNORECASE)
_BLOCK_WIDTH: Final = 60
_CODON_PEP_WIDTH: Final = 20
_MIN_ID_WIDTH: Final = 10


def chunk(seq: str, width: int) -> list[str]:
    """pal2nal.pl sub my1while (line 1920) -- fixed-width chunking."""
    return [seq[i : i + width] for i in range(0, len(seq), width)] or [""]


@dataclass
class Alignment:
    ids: list[str]
    codonaln: list[str] = field(default_factory=list)
    coloraln: list[str] = field(default_factory=list)
    maskseq: str = ""


def build(
    ids: list[str],
    aaseqs: list[str],
    codonseqs: list[str],
    mismatches: set[tuple[str, int]],
    blockseq: str,
    *,
    blockonly: bool = False,
    nomismatch: bool = False,
) -> Alignment:
    """pal2nal.pl:736-838 -- reinsert gaps and apply the column filters."""
    aln = Alignment(ids=ids, codonaln=["" for _ in ids], coloraln=["" for _ in ids])
    cursors = [0 for _ in ids]
    has_blocks = "#" in blockseq
    alnlen = len(aaseqs[0]) if aaseqs else 0

    for col in range(alnlen):
        column = [seq[col] if col < len(seq) else "-" for seq in aaseqs]

        # the column is as wide as its widest entry needs
        width = 3
        for aa in column:
            if _DIGIT.fullmatch(aa):
                width = max(width, -(-int(aa) // 3) * 3)
        put = True
        if has_blocks and blockonly and blockseq[col : col + 1] != "#":
            put = False
        if nomismatch and any((sid, col) in mismatches for sid in ids):
            put = False

        for k, aa in enumerate(column):
            if _DIGIT.fullmatch(aa):
                take = int(aa)
                piece = codonseqs[k][cursors[k] : cursors[k] + take]
                cursors[k] += take
                piece += "-" * (width - len(piece))
                colour = "-" * width
            elif _GAP.fullmatch(aa):
                piece = "-" * width
                colour = "-" * width
            else:
                piece = codonseqs[k][cursors[k] : cursors[k] + 3]
                cursors[k] += 3
                piece += "-" * (width - len(piece))
                mark = "R" if (ids[k], col) in mismatches else "-"
                colour = mark * 3 + "-" * (width - 3)
            if put:
                aln.codonaln[k] += piece
                aln.coloraln[k] += colour
        if put and not blockonly:
            marker = blockseq[col : col + 1] or " "
            aln.maskseq += marker * width
    return aln


def remove_gaps(aln: Alignment) -> Alignment:
    """pal2nal.pl:845-873 -- drop any codon column that is gapped in any
    sequence or that holds an in-frame stop."""
    keep: list[int] = []
    length = len(aln.codonaln[0]) if aln.codonaln else 0
    for pos in range(0, length - 2, 3):
        ok = True
        for seq in aln.codonaln:
            codon = seq[pos : pos + 3]
            if "-" in codon or _STOP.fullmatch(codon):
                ok = False
                break
        if ok:
            keep.append(pos)
    out = Alignment(ids=aln.ids)
    out.codonaln = ["".join(s[p : p + 3] for p in keep) for s in aln.codonaln]
    out.coloraln = ["".join(s[p : p + 3] for p in keep) for s in aln.coloraln]
    out.maskseq = "".join(aln.maskseq[p : p + 3] for p in keep)
    return out


def _emit(out: TextIO, text: str, colour: str, html: bool) -> None:
    """One chunk, red-flagging mismatched characters individually in html
    mode exactly as v14 does (one FONT element per character)."""
    if html and "R" in colour:
        for i, ch in enumerate(text):
            if colour[i : i + 1] == "R":
                out.write(f"<FONT color='red'>{ch}</FONT>")
            else:
                out.write(ch)
        out.write("\n")
    else:
        out.write(text + "\n")


def write_alignment(
    out: TextIO,
    aln: Alignment,
    aaseqs: list[str],
    *,
    outform: str,
    html: bool,
    blockonly: bool,
    blockseq: str,
) -> None:
    idw = max([_MIN_ID_WIDTH] + [len(i) for i in aln.ids])
    show_mask = not blockonly and "#" in blockseq

    if outform == "paml":
        out.write(f" {len(aln.ids):3d} {len(aln.codonaln[0]):6d}\n")
        for k, sid in enumerate(aln.ids):
            out.write(sid + "\n")
            for text, colour in zip(
                chunk(aln.codonaln[k], _BLOCK_WIDTH), chunk(aln.coloraln[k], _BLOCK_WIDTH)
            ):
                _emit(out, text, colour, html)
        return

    if outform == "fasta":
        for k, sid in enumerate(aln.ids):
            out.write(">" + sid + "\n")
            for text, colour in zip(
                chunk(aln.codonaln[k], _BLOCK_WIDTH), chunk(aln.coloraln[k], _BLOCK_WIDTH)
            ):
                _emit(out, text, colour, html)
        return

    if outform == "codon":
        _write_codon(out, aln, aaseqs, idw, html, show_mask)
        return

    out.write("CLUSTAL W multiple sequence alignment\n")
    out.write("\n")
    blocks = chunk(aln.codonaln[0], _BLOCK_WIDTH)
    masks = chunk(aln.maskseq, _BLOCK_WIDTH)
    for b in range(len(blocks)):
        for k, sid in enumerate(aln.ids):
            out.write(f"{sid:<{idw}}    ")
            _emit(
                out,
                chunk(aln.codonaln[k], _BLOCK_WIDTH)[b],
                chunk(aln.coloraln[k], _BLOCK_WIDTH)[b],
                html,
            )
        if show_mask:
            mask = masks[b] if b < len(masks) else ""
            out.write(f"{' ':<{idw}}    {mask}\n")
        out.write("\n")


def _peptide_rows(aaseqs: list[str]) -> list[str]:
    """pal2nal.pl:959-991 -- one display slot per alignment column, widened
    where a frame shift consumes more than three nucleotides.

    v14 sized the widened column from a loop variable left over from the
    scan for the column maximum, i.e. from whichever residue the *last*
    sequence happened to have there. When the frame-shift digit belonged to
    any other sequence the peptide row came out a column short and drifted
    against the codon row. The column maximum is used here instead.
    """
    if not any(_DIGIT.search(s) for s in aaseqs):
        return list(aaseqs)
    rows = ["" for _ in aaseqs]
    alnlen = len(aaseqs[0]) if aaseqs else 0
    for col in range(alnlen):
        column = [s[col] if col < len(s) else "-" for s in aaseqs]
        maxaan = 0
        for aa in column:
            if _DIGIT.fullmatch(aa):
                maxaan = max(maxaan, int(aa))
        slots = (maxaan - 1) // 3 + 1 if maxaan >= 4 else 1
        for k, aa in enumerate(column):
            rows[k] += aa + "-" * (slots - 1)
    return rows


def _write_codon(
    out: TextIO, aln: Alignment, aaseqs: list[str], idw: int, html: bool, show_mask: bool
) -> None:
    peps = _peptide_rows(aaseqs)
    blocks = chunk(aln.codonaln[0], _BLOCK_WIDTH)
    masks = chunk(aln.maskseq, _BLOCK_WIDTH)
    for b in range(len(blocks)):
        for k, sid in enumerate(aln.ids):
            pep_chunks = chunk(peps[k], _CODON_PEP_WIDTH)
            pep = pep_chunks[b] if b < len(pep_chunks) else ""
            out.write(f"{'':<{idw}}    ")
            out.write("   ".join(pep) + "\n")

            codon_chunk = chunk(aln.codonaln[k], _BLOCK_WIDTH)[b]
            colour_chunk = chunk(aln.coloraln[k], _BLOCK_WIDTH)[b]
            out.write(f"{sid:<{idw}}    ")
            spaced = " ".join(chunk(codon_chunk, 3))
            spaced_colour = " ".join(chunk(colour_chunk, 3))
            _emit(out, spaced, spaced_colour, html)
        if show_mask:
            mask = masks[b] if b < len(masks) else ""
            out.write(f"{' ':<{idw}}    {' '.join(chunk(mask, 3))}\n")
        out.write("\n")
