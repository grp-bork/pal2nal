# PAL2NAL

Convert a protein sequence alignment and the corresponding DNA (or mRNA)
sequences into a codon-based DNA alignment.

The codon assignment tolerates mismatches between the protein and the DNA,
untranslated regions, polyA tails, in-frame stop codons and frame shifts,
which makes it usable on pseudogenes as well as on genes.

This is **v16**, a Python port of `pal2nal.pl` v14 by Mikita Suyama. It
reproduces v14 byte for byte except where [CHANGELOG.md](CHANGELOG.md)
records a change — most visibly [`-partial`](#partial), which converts a
peptide and a CDS that do not fully correspond instead of aborting the run.
Conversion is also faster than v14 — roughly 3× on a large alignment and
10–14× on the inputs that made v14 slowest — because the codon
patterns are matched directly instead of being compiled into a peptide-length
regular expression for every sequence. On a small alignment the two are
comparable: Python's ~50 ms of startup outweighs the work.

> Suyama M, Torrents D, Bork P (2006). PAL2NAL: robust conversion of
> protein sequence alignments into the corresponding codon alignments.
> *Nucleic Acids Research* **34**: W609–W612.

## Install

```sh
pip install git+https://github.com/grp-bork/pal2nal
```

Python 3.11 or newer. No runtime dependencies.

## Use

```sh
pal2nal pep.aln nuc.fasta [nuc.fasta ...] [options] > out.codon
```

Try it on the bundled example:

```sh
pal2nal examples/test.aln examples/test.nuc
```

`pep.aln` is a protein alignment in CLUSTAL or FASTA format; the format is
detected automatically. The DNA may be one multi-FASTA file or several.
Sequences are paired by ID when the IDs correspond and by input order
otherwise.

| Option | Effect |
|---|---|
| `-output clustal\|paml\|fasta\|codon` | output format, default `clustal` |
| `-nogap` | drop columns with gaps and in-frame stop codons |
| `-nomismatch` | drop codons that disagree with the protein |
| `-partial` | convert what can be matched when a peptide and its DNA disagree, instead of aborting |
| `-blockonly` | keep only columns marked `#` under the alignment |
| `-codontable N` | NCBI genetic code, default 1 |
| `-html` | HTML output |
| `-nostderr` | suppress messages |
| `-h` | usage |
| `--version` | print the version and stop |

Supported codon tables: 1–6, 9–16, 21–33 (every genetic code NCBI
defines), plus 34–37. NCBI assigns no ID past 33; 34–37 are the alternative
bacterial codes of Shulgina & Eddy (2021) under Wikipedia's numbering, which
a future NCBI assignment could collide with. Tables 7 and 8 were retired by
NCBI in 1995 and 17–20 have never been assigned, so neither is accepted.

Mark in-frame stop codons with `*` or `_`, and frame shifts with a digit
giving the number of nucleotides the column consumes.

<a id="partial"></a>

## `-partial`

Without it, a single peptide whose DNA does not correspond to it aborts the
whole run — every other sequence in the alignment is discarded along with the
offending one:

```
#---  ERROR: inconsistency between the following pep and nuc seqs  ---#
```

Four common shapes of real Ensembl and NCBI data reach that error: an intron
left in the DNA, an indel of any size between the peptide and its CDS, and a
CDS truncated at either end. `-partial` places each ten-residue anchor against
the DNA independently, keeps the codons it can, gaps the residues it cannot,
and carries on. An intron is recovered exactly — the codons either side of it
are the true CDS, byte for byte.

The flag is opt-in and engages only after the normal matching has failed, so
it changes nothing about an input that already converts. When it does engage
it reports, per sequence, how much was placed:

```
#  PARTIAL: seq1 60/60 residues placed in 2 segments
#  PARTIAL: seq2 50/60 residues placed
#  UNMATCHED: seq3 0/60 residues placed
```

The exit status stays 0 — not aborting is the point — so that report is
how a caller learns what happened. `-nostderr` silences it; `-nomismatch` and
`-blockonly` do not, since those choose which codons are shown and must not
hide how much of a sequence was matched. A residue with no codon counts as a
mismatched column: `-nomismatch` drops it, `-html` marks it, `-nogap` removes
it as the gap it is.

It is deliberately not a splice-aware aligner. Anchors are chained greedily
from left to right, so a tandem repeat can attract one to the wrong copy; the
residue counts in the report are what shows that. Validation is untouched —
a ragged alignment, a duplicate ID or a bad alphabet is still refused.

## Tests

```sh
python -m venv .venv && .venv/bin/pip install -e '.[dev]'
.venv/bin/python -m pytest
```

The suite compares the port against output captured from the original Perl.
`tests/golden/` holds what `pal2nal.v14.pl` produces for every case in
`tests/cases.tsv`; `tests/expected/` holds the intended output for the cases
where the port deliberately differs — 75 of the 171 cases, each registered in
`tests/divergences.tsv` with a reason. **Those files are evidence: never edit
one to make a test pass.** See [tests/README.md](tests/README.md).

`tests/reference/pal2nal.v14.pl` is the original implementation, kept so the
corpus can be regenerated (`cd tests && sh generate_golden.sh`, needs `perl`)
and for reproducibility of previous results.

## Repository layout

| Path | Contents |
|---|---|
| `pal2nal/` | the Python package |
| `tests/` | the golden corpus and the test suite |
| `examples/` | a small alignment to try the tool on |
| `CHANGELOG.md` | what changed since v14 and why |
| `PORTING.md` | the porting decisions and the evidence behind them |

The web front end, the full Perl and CGI history of the original EMBL
service, and the notes on how that history was reconstructed live in a
separate archive repository and are deliberately not part of this package.

## Licence

MIT — see [LICENSE](LICENSE).
