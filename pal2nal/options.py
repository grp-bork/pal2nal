"""Command-line option parsing, faithful to pal2nal.v14.pl lines 245-324.

The original parses argv in a single left-to-right pass through an if/elsif
chain. Two options take a value, and they take it by setting a flag that
makes the *next* token be consumed as the value whatever it looks like --
so ``-output -html`` reads "-html" as the format and fails. A trailing
``-output`` with nothing after it is silently ignored. Both behaviours are
reproduced here; the golden corpus pins them.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from typing import Final, TextIO

from .codontables import SUPPORTED

OUTPUT_FORMATS: Final[tuple[str, ...]] = ("clustal", "paml", "fasta", "codon")

_LEADING_NUMBER = re.compile(r"\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)")


class Exit(Exception):
    """Raised where the Perl calls exit, so the CLI stays testable."""

    def __init__(self, status: int = 0) -> None:
        super().__init__(status)
        self.status = status


@dataclass
class Options:
    alnfile: str | None = None
    nucfiles: list[str] = field(default_factory=list)
    outform: str = "clustal"
    nogap: bool = False
    nomismatch: bool = False
    blockonly: bool = False
    html: bool = False
    nostderr: bool = False
    codontable: int = 1


def perl_num(text: str) -> float:
    """Perl's loose string-to-number coercion, as used by the `!=` tests
    that validate -codontable: " 1", "1.0" and "1abc" all read as 1, and
    anything with no leading numeric part reads as 0."""
    m = _LEADING_NUMBER.match(text)
    return float(m.group(1)) if m else 0.0


def parse_args(argv: list[str], *, stderr: TextIO | None = None) -> Options:
    """Parse argv. Raises Exit wherever the Perl exits."""
    err = stderr if stderr is not None else sys.stderr

    # pal2nal.pl:245 -- fewer than two arguments prints usage and stops
    if len(argv) < 2:
        from .help import showhelp

        showhelp(err)
        raise Exit(0)

    opt = Options()
    get_outform = False
    get_codontable = False
    valid_tables = {float(n) for n in SUPPORTED}

    for arg in argv:
        if arg == "-h":
            # v14 printed usage here and carried on converting; v15 stops
            from .help import showhelp

            showhelp(err)
            raise Exit(0)
        elif arg == "-output":
            get_outform = True
        elif get_outform:
            if arg.startswith("-"):
                err.write("\nERROR:  -output needs a format: clustal, paml, fasta, or codon\n\n")
                raise Exit(1)
            opt.outform = arg
            if opt.outform not in OUTPUT_FORMATS:
                err.write("\nERROR:  valid output format: clustal, paml, fasta, or codon\n\n")
                raise Exit(1)
            get_outform = False
        elif arg == "-blockonly":
            opt.blockonly = True
        elif arg == "-nogap":
            opt.nogap = True
        elif arg == "-nomismatch":
            opt.nomismatch = True
        elif arg == "-codontable":
            get_codontable = True
        elif get_codontable:
            if arg.startswith("-"):
                err.write("\nERROR:  -codontable needs a table number\n\n")
                raise Exit(1)
            if perl_num(arg) not in valid_tables:
                err.write(f"\nERROR:  invalid codontable number, {arg}!!\n\n")
                raise Exit(1)
            opt.codontable = int(perl_num(arg))
            get_codontable = False
        elif arg == "-html":
            opt.html = True
        elif arg == "-nostderr":
            opt.nostderr = True
        elif not opt.alnfile:
            # truthiness, as in the Perl: a file literally named "0" is
            # falsy and would be overwritten by the next positional
            opt.alnfile = arg
        else:
            opt.nucfiles.append(arg)

    if get_outform:
        err.write("\nERROR:  -output needs a format: clustal, paml, fasta, or codon\n\n")
        raise Exit(1)
    if get_codontable:
        err.write("\nERROR:  -codontable needs a table number\n\n")
        raise Exit(1)

    return opt


def check_combinations(opt: Options, out: TextIO, err: TextIO) -> None:
    """pal2nal.pl:316 -- "codon" output is incompatible with the filters.

    Runs after the <pre> banner has already been emitted in html mode, and
    reports on stdout there rather than stderr.
    """
    if opt.outform == "codon" and (opt.blockonly or opt.nogap or opt.nomismatch):
        if opt.html:
            out.write(
                '\nERROR:  if "codon(Output format)" is selected, '
                'don\'t use "Remove gaps, inframe stop codons" or '
                '"Remove mismatches" or "Use only selected positions".\n\n'
            )
        else:
            err.write(
                '\nERROR:  "-outform codon" is not valid with '
                "-blockonly, -nogap, -nomismatch\n\n"
            )
        raise Exit(1)
