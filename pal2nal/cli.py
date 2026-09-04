"""Command-line entry point: the main body of pal2nal.v14.pl.

Stream routing follows v14: warnings go to stderr, or to stdout wrapped in
<pre> under -html, and are suppressed entirely by -nostderr or -nomismatch.
"""

from __future__ import annotations

import html as html_module
import sys
from typing import Final, TextIO

from . import __version__, inputs, output
from .convert import NO_MATCH, pn2codon
from .options import Exit, Options, check_combinations, parse_args
from .validate import InputError

RULE = "#------------------------------------------------------------------------#"

#: New in v15; v14 had no way to report its own version. Both spellings are
#: accepted: the tool's own options take a single dash, but "--version" is
#: what everything else in a pipeline tries first.
_VERSION_FLAGS: Final = frozenset({"-version", "--version"})


def _esc(text: str, html: bool) -> str:
    """v14 wrote ids and sequences into -html output unescaped, so an id
    containing "<" or "&" produced broken markup. v15 escapes them."""
    return html_module.escape(text, quote=False) if html else text


def _fail(opt: Options, out: TextIO, err: TextIO, message: str) -> None:
    """Report an error from inside run() and stop.

    The stream is v14's: stderr, or stdout under -html. Under -html the
    <pre> opened at the top of run() is closed first, so a rejected input
    cannot leave the page with unbalanced markup swallowing everything
    after it.
    """
    if opt.html:
        out.write(message)
        out.write("</pre>\n")
    else:
        err.write(message)
    raise Exit(1)


def run(opt: Options, out: TextIO, err: TextIO) -> int:
    if opt.html:
        out.write("<pre>\n")
    check_combinations(opt, out, err)

    try:
        nuc = inputs.read_nucleotides(opt.nucfiles)
        aln = inputs.read_alignment(opt.alnfile) if opt.alnfile else inputs.Alignment()
    except InputError as exc:
        # the message carries no input text, but a sequence id may appear
        # in it, so it is escaped like everything else under -html
        _fail(opt, out, err, _esc(exc.message, opt.html))

    # v15: a repeated FASTA header used to be appended to the id list twice
    # while both records' residues merged under one key
    duplicates = [i for i in set(aln.ids) if aln.ids.count(i) > 1]
    if duplicates:
        _fail(
            opt, out, err,
            f"\nERROR: duplicate sequence ID(s) in {opt.alnfile}: "
            f"{', '.join(sorted(duplicates))}\n\n",
        )

    if len(aln.ids) != len(nuc.ids):
        naa, nnuc = len(aln.ids), len(nuc.ids)
        msg = f"\nERROR: number of input seqs differ (aa: {naa};  nuc: {nnuc})!!\n\n"
        if not opt.html:
            # v14 lists the ids only on stderr
            msg += "   aa  '{}'\n".format(" ".join(aln.ids))
            msg += "   nuc '{}'\n".format(" ".join(nuc.ids))
        _fail(opt, out, err, msg)

    aaseqs = [inputs.apply_frameshift_markers(s) for s in aln.sequences]
    correspondence = inputs.id_correspondence(aln.ids, nuc.ids)

    codonseqs: list[str] = []
    messages: list[str] = []
    unknown: list[str] = []
    mismatches: set[tuple[str, int]] = set()

    for i, aa_id in enumerate(aln.ids):
        nuc_id = aa_id if correspondence == "sameID" else nuc.ids[i]
        res = pn2codon(aaseqs[i], nuc.id2seq.get(nuc_id, ""), opt.codontable)
        if res.result == NO_MATCH:
            _fail(
                opt, out, err,
                _inconsistency(opt, aa_id, nuc_id, aaseqs[i], nuc.id2seq.get(nuc_id, "")),
            )
        codonseqs.append(res.codonseq)
        for message in res.unknown:
            unknown.append(f"WARNING: {aa_id} {message}")
        for message in res.messages:
            messages.append(f"WARNING: {aa_id} {message}")
        # a column is a mismatch either way: an unknown character was taken
        # as X, so its codon was never verified. The raw text is kept for
        # position parsing here; escaping happens at write time, so -html
        # output cannot carry markup smuggled in through an id or a residue
        for message in res.unknown + res.messages:
            parts = message.split()
            if len(parts) > 1 and parts[1].rstrip(":").isdigit():
                mismatches.add((aa_id, int(parts[1].rstrip(":")) - 1))

    if aln.blockseq and "#" in aln.blockseq and opt.blockonly:
        kept = []
        for message in messages:
            parts = message.split()
            pos = parts[3].rstrip(":") if len(parts) > 3 else ""
            if pos.isdigit() and aln.blockseq[int(pos) - 1 : int(pos)] == "#":
                kept.append(message)
        messages = kept

    # v15: -nomismatch and -blockonly select which *codons* are reported,
    # so neither may hide a character that is not a residue at all. Only
    # -nostderr, which asks for silence outright, still suppresses those.
    reported = list(unknown)
    if not opt.nomismatch:
        if opt.codontable != 1:
            reported.insert(0, f"Codontable {opt.codontable} is used")
        reported += messages

    if not opt.nostderr and reported:
        sink = out if opt.html else err
        sink.write(RULE + "\n")
        if not opt.html:
            sink.write("#  Input files:  {} {}\n".format(opt.alnfile, " ".join(opt.nucfiles)))
        for message in reported:
            sink.write(f"#  {_esc(message, opt.html)}\n")
        sink.write(RULE + "\n\n")

    built = output.build(
        aln.ids, aaseqs, codonseqs, mismatches, aln.blockseq,
        blockonly=opt.blockonly, nomismatch=opt.nomismatch,
    )
    if opt.nogap:
        built = output.remove_gaps(built, opt.codontable)
    built.ids = [_esc(i, opt.html) for i in built.ids]

    output.write_alignment(
        out, built, aaseqs,
        outform=opt.outform, html=opt.html,
        blockonly=opt.blockonly, blockseq=aln.blockseq,
    )
    if opt.html:
        out.write("</pre>\n")
    return 0


def _inconsistency(opt: Options, aa_id: str, nuc_id: str, pep: str, dna: str) -> str:
    """pal2nal.pl:538-660, without the bl2seq diagnostic (dropped: the tool
    is retired and the message it produced was only advisory).

    This is the one message that does repeat the input back. It is reached
    only after validate.py has passed both sequences, so what it echoes is
    printable ASCII, and it is escaped under -html.
    """
    parts = ["#---  ERROR: inconsistency between the following pep and nuc seqs  ---#\n"]
    parts.append(f">{_esc(aa_id, opt.html)}\n")
    parts += [_esc(line, opt.html) + "\n" for line in output.chunk(pep.replace("-", ""), 60)]
    parts.append(f">{_esc(nuc_id, opt.html)}\n")
    parts += [_esc(line, opt.html) + "\n" for line in output.chunk(dna, 60)]
    parts.append(RULE + "\n\n")
    return "".join(parts)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    out, err = sys.stdout, sys.stderr

    # Before anything else, including the "fewer than two arguments prints
    # usage" rule -- "pal2nal --version" is a single argument and has to
    # report the version rather than the usage text. On stdout, unlike -h:
    # this exists to be read by whatever is recording which version ran.
    if _VERSION_FLAGS & set(args):
        out.write(f"pal2nal {__version__}\n")
        return 0

    try:
        opt = parse_args(args, stderr=err)
        if not opt.alnfile or not opt.nucfiles:
            err.write("\nERROR: an alignment file and at least one DNA file are required\n\n")
            return 1
        return run(opt, out, err)
    except Exit as exc:
        return exc.status
    except OSError as exc:
        err.write("Can't open %s\n" % (exc.filename or ""))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
