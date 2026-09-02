"""Input validation: what the converter refuses to look at.

v14 read whatever bytes it was given. Anything that was not a residue
became "X" on the peptide side, and on the DNA side anything that was not
a letter was silently deleted while a stray letter was kept and quietly
broke the codon match -- the run then died with the generic "inconsistency
between the following pep and nuc seqs", which dumps both sequences back
at the caller. Under -html that dump is the web page.

Three gates are applied before any of that:

* non-ASCII input is refused, and the message never repeats the offending
  text -- only how much there is and where the first one is. Bytes that
  are not valid UTF-8 arrive as surrogates (inputs.py reads with
  errors="surrogateescape") and are non-ASCII too, so a binary file is
  rejected here rather than being pattern-matched against;
* ASCII control characters other than tab, CR and LF are refused the same
  way, named by code point rather than written out, so nothing invisible
  or terminal-controlling reaches stdout or an HTML page;
* DNA characters outside the IUPAC alphabet are refused by name. They are
  printable ASCII by the time this runs, so quoting them is safe, and
  saying "'E' at nucleotide 42" beats the sequence dump it replaces.

Peptide residues are deliberately *not* gated here: an unknown residue has
a defined fallback (it is taken as X) and convert.py already reports each
one with its alignment position, exactly as v14 did. An unknown nucleotide
has no fallback, so it is fatal.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Final

#: the same three line endings inputs.normalise_newlines collapses; kept
#: separate so this module has no import back into inputs
_NEWLINES: Final = re.compile(r"\x0D\x0A|\x0D|\x0A")

_NON_ASCII: Final = re.compile(r"[^\x00-\x7F]")
_CONTROL: Final = re.compile(r"[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]")

#: the four bases, U for RNA input, and the IUPAC ambiguity codes
NUCLEOTIDES: Final = "ACGTU"
AMBIGUITY: Final = "RYSWKMBDHVN"
_NUC_ALPHABET: Final = frozenset(NUCLEOTIDES + AMBIGUITY)

#: at most this many sequences are named in one alphabet complaint, so a
#: thoroughly corrupt file produces a report and not a second data dump
_MAX_REPORTED: Final = 5


class InputError(Exception):
    """A rejected input.

    ``message`` is ready to print and contains no input text beyond, at
    most, a sequence id and individual offending characters that are known
    to be printable ASCII.
    """

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


def _position(text: str, index: int) -> tuple[int, int]:
    """1-based line and column of `index`, counting CRLF, CR and LF."""
    prefix = text[:index]
    starts = [m.end() for m in _NEWLINES.finditer(prefix)]
    return len(starts) + 1, index - (starts[-1] if starts else 0) + 1


def check_text(text: str, source: str) -> None:
    """Refuse non-ASCII and control characters in one raw input file.

    `source` names the input in the message -- a path from the CLI, a form
    field from the web front end, which must not leak its scratch paths.
    """
    match = _NON_ASCII.search(text)
    if match:
        count = len(_NON_ASCII.findall(text))
        line, column = _position(text, match.start())
        raise InputError(
            f"\nERROR: {source} contains {count} non-ASCII character"
            f"{'' if count == 1 else 's'} "
            f"(first at line {line}, column {column})\n"
            "       expected plain ASCII; the offending text is not echoed\n\n"
        )

    match = _CONTROL.search(text)
    if match:
        count = len(_CONTROL.findall(text))
        line, column = _position(text, match.start())
        raise InputError(
            f"\nERROR: {source} contains {count} control character"
            f"{'' if count == 1 else 's'} "
            f"(first 0x{ord(match.group(0)):02X} at line {line}, column {column})\n"
            "       only tab, carriage return and newline are allowed\n\n"
        )


def _offenders(seq: str) -> list[tuple[str, int, int]]:
    """Characters of `seq` outside the nucleotide alphabet, in order of
    first appearance, as (character, count, first 1-based position)."""
    found: dict[str, list[int]] = {}
    for i, ch in enumerate(seq):
        if ch.upper() not in _NUC_ALPHABET:
            found.setdefault(ch.upper(), [0, i + 1])[0] += 1
    return [(ch, n, first) for ch, (n, first) in found.items()]


def check_nucleotide_alphabet(id2seq: Mapping[str, str]) -> None:
    """Refuse DNA sequences holding characters that are not nucleotides.

    Only letters can reach this: inputs.read_nucleotides strips everything
    else, because numbered and gapped FASTA variants rely on that.
    """
    bad = [(i, _offenders(s)) for i, s in id2seq.items()]
    bad = [(i, o) for i, o in bad if o]
    if not bad:
        return

    lines = [
        "\nERROR: unexpected characters in the DNA sequence(s):\n",
    ]
    for seq_id, offenders in bad[:_MAX_REPORTED]:
        detail = ", ".join(
            f"'{ch}' (x{n}, first at nucleotide {first})" for ch, n, first in offenders
        )
        lines.append(f"       {seq_id}: {detail}\n")
    if len(bad) > _MAX_REPORTED:
        lines.append(f"       ... and {len(bad) - _MAX_REPORTED} more sequence(s)\n")
    lines.append(
        f"       expected {', '.join(NUCLEOTIDES)} or an IUPAC ambiguity code "
        f"({' '.join(AMBIGUITY)})\n\n"
    )
    raise InputError("".join(lines))
