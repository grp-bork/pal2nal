"""Reading the peptide alignment and the DNA sequences.

Mirrors pal2nal.v14.pl lines 330-453. Two details here are accidents of
Perl that the output depends on, so they are reproduced deliberately:

* The alignment format is sniffed from the file read as a *single string*.
  The original does `undef $/` at line 242, which turns the `while
  (<ALNFILE>)` at line 366 into one iteration holding the whole file, and
  the `^` anchors then only match at the very start of it. A file whose
  first byte is a newline therefore fails every test and falls through to
  the "clustal" default -- a FASTA alignment with a leading blank line is
  parsed as CLUSTAL and comes apart. Verified against the Perl.

* The CLUSTAL conservation/marker line is identified by *position*, never
  by content: whatever line follows a sequence line is treated as the
  marker line for that block, sliced at the same columns.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .validate import check_alignment_widths, check_nucleotide_alphabet, check_text

_NEWLINES = re.compile(r"\x0D\x0A|\x0D|\x0A")
_FASTA_HEADER = re.compile(r"^>(\S+)")
_ID_THEN_SPACE = re.compile(r"^\S+\s+")
_NON_LETTER = re.compile(r"[^a-zA-Z]")
_WHITESPACE = re.compile(r"\s+")


def normalise_newlines(data: str) -> str:
    """CRLF, bare CR and LF all become LF (pal2nal.pl:338, v11/v13)."""
    return _NEWLINES.sub("\n", data)


def perl_split_lines(data: str) -> list[str]:
    """Perl's split(/\\n/, $data), which drops trailing empty fields."""
    parts = data.split("\n")
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def detect_alignment_type(data: str) -> str | None:
    """Decide the alignment format from the first meaningful line.

    v14 tested the file read as one string, so the anchors only matched at
    byte 0 and a leading blank line sent a FASTA alignment down the CLUSTAL
    path, where it came apart. v15 skips blank and "#" lines first.

    Returns None for a file with no meaningful line at all.
    """
    for line in data.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith("CLUSTAL"):
            return "clustal"
        if line.startswith(">"):
            return "fasta"
        if line.startswith("Gblocks"):
            return "gblocks"
        return "clustal"
    return None


@dataclass
class Alignment:
    ids: list[str] = field(default_factory=list)
    id2aln: dict[str, str] = field(default_factory=dict)
    #: the marker-line characters under the alignment, used by -blockonly
    blockseq: str = ""

    @property
    def sequences(self) -> list[str]:
        """pal2nal.pl:441 -- looked up by id, so a duplicated FASTA header
        yields the same merged sequence twice."""
        return [self.id2aln.get(i, "") for i in self.ids]


def read_alignment(path: str) -> Alignment:
    with open(path, encoding="utf-8", errors="surrogateescape") as fh:
        data = fh.read()
    # undecodable bytes are surrogates here, so this rejects them too
    check_text(data, path)

    aln_type = detect_alignment_type(data)
    aln = Alignment()
    if aln_type is None:
        return aln

    data = normalise_newlines(data)
    lines = perl_split_lines(data)

    if aln_type == "fasta":
        _parse_fasta(lines, aln)
    elif aln_type == "gblocks":
        _parse_gblocks(lines, aln)
    else:
        _parse_clustal(lines, aln)
    check_alignment_widths(aln.ids, aln.sequences)
    return aln


def _parse_clustal(lines: list[str], aln: Alignment) -> None:
    """pal2nal.pl:390-411."""
    seen: dict[str, int] = {}
    getblock = False
    tmplen = idspc = subalnlen = 0

    for line in lines:
        is_data = (
            not line.startswith("CLUSTAL")
            and bool(line[:1].strip())      # /^\S+/
            and not line.startswith("#")
        )
        if is_data:
            line = re.sub(r"\s+$", "", line)
            dat = _WHITESPACE.split(line)
            seq_id = dat[0]
            seen[seq_id] = seen.get(seq_id, 0) + 1
            if seen[seq_id] == 1:
                aln.ids.append(seq_id)
            chunk = dat[1].upper() if len(dat) > 1 else ""
            aln.id2aln[seq_id] = aln.id2aln.get(seq_id, "") + chunk

            tmplen = len(line)
            m = _ID_THEN_SPACE.match(line)
            idspc = len(m.group(0)) if m else 0
            subalnlen = len(chunk)
            getblock = True
        elif getblock:
            # the marker line, taken purely by position
            padded = line + " " * max(0, tmplen - len(line))
            aln.blockseq += padded[idspc : idspc + subalnlen]
            getblock = False


def _parse_fasta(lines: list[str], aln: Alignment) -> None:
    """pal2nal.pl:412-420. Note the missing dedup: a repeated header is
    appended to ids again, while the sequence merges under one key."""
    tmpid = ""
    for line in lines:
        m = _FASTA_HEADER.match(line)
        if m:
            tmpid = m.group(1)
            aln.ids.append(tmpid)
        else:
            seq = _WHITESPACE.sub("", line).upper()
            aln.id2aln[tmpid] = aln.id2aln.get(tmpid, "") + seq


def _parse_gblocks(lines: list[str], aln: Alignment) -> None:
    """pal2nal.pl:421-438. Reachable only for a file literally starting
    with "Gblocks"; the -a option that selected it went away in v8."""
    seen: dict[str, int] = {}
    getaln = False
    for line in lines:
        if re.match(r"^\s+=", line):
            getaln = True
            continue
        if line.startswith("Parameters"):
            getaln = False
            continue
        if not getaln or not line.strip():
            continue
        dat = _WHITESPACE.split(line.rstrip())
        if line.startswith("Gblocks"):
            if len(dat) > 1:
                aln.blockseq += dat[1]
        else:
            seq_id = dat[0]
            seen[seq_id] = seen.get(seq_id, 0) + 1
            if seen[seq_id] == 1:
                aln.ids.append(seq_id)
            if len(dat) > 1:
                aln.id2aln[seq_id] = aln.id2aln.get(seq_id, "") + dat[1].upper()


def apply_frameshift_markers(seq: str) -> str:
    """pal2nal.pl:450-453 -- the tfastx/tfasty frame-shift notations.

    A backslash becomes "1". A slash, any following gaps, and the residue
    after them become those gaps plus "2": the residue is consumed.
    """
    seq = seq.replace("\\", "1")
    return re.sub(r"/(-*)[A-Z*]", r"-\g<1>2", seq)


@dataclass
class NucSequences:
    ids: list[str] = field(default_factory=list)
    id2seq: dict[str, str] = field(default_factory=dict)


def read_nucleotides(paths: list[str]) -> NucSequences:
    """pal2nal.pl:330-351. Accumulates across all files without resetting,
    so an id repeated in another file simply extends the same sequence.
    Case is preserved here; only the peptide side is upper-cased."""
    nuc = NucSequences()
    # accumulated per id and joined once: "seq = seq + line" over a FASTA
    # file is quadratic in the length of the sequence
    chunks: dict[str, list[str]] = {}
    tmpid = ""
    for path in paths:
        with open(path, encoding="utf-8", errors="surrogateescape") as fh:
            data = fh.read()
        check_text(data, path)
        for line in perl_split_lines(normalise_newlines(data)):
            if line.startswith("#") or not line.strip():
                continue
            m = _FASTA_HEADER.match(line)
            if m:
                tmpid = m.group(1)
                nuc.ids.append(tmpid)
            else:
                chunks.setdefault(tmpid, []).append(line)
    # stripping per sequence rather than per line: the substitution is
    # per character either way, and there is one call instead of thousands
    nuc.id2seq = {
        seq_id: _NON_LETTER.sub("", "".join(parts)) for seq_id, parts in chunks.items()
    }
    # only letters survive the strip above, so this catches stray residues
    # and prose, not the digits and dashes of numbered or gapped FASTA
    check_nucleotide_alphabet(nuc.id2seq)
    return nuc


def common_elements(first: list[str], second: list[str]) -> list[str]:
    """pal2nal.pl:480 sub common_elem -- entries of `second` that occur
    anywhere in `first`, keeping `second`'s order and its duplicates."""
    mark = set(first)
    return [x for x in second if x in mark]


def id_correspondence(aa_ids: list[str], nuc_ids: list[str]) -> str:
    """Use id-based pairing only for a genuine one-to-one match.

    v14 compared list lengths (`$#commonids == $#aaid`), which duplicates
    can satisfy without the id sets matching, so it could pick id-based
    pairing and then pair a peptide with a missing DNA sequence.
    """
    if len(aa_ids) != len(nuc_ids):
        return "ordered"
    if len(set(aa_ids)) != len(aa_ids) or len(set(nuc_ids)) != len(nuc_ids):
        return "ordered"
    return "sameID" if set(aa_ids) == set(nuc_ids) else "ordered"
