# PAL2NAL

Convert a protein sequence alignment and the corresponding DNA (or mRNA)
sequences into a codon-based DNA alignment.

The codon assignment tolerates mismatches between the protein and the DNA,
untranslated regions, polyA tails, in-frame stop codons and frame shifts,
which makes it usable on pseudogenes as well as on genes.

This is **v15**, a Python port of `pal2nal.pl` v14 by Mikita Suyama. It
reproduces v14 byte for byte except where [CHANGELOG.md](CHANGELOG.md)
records a change. Conversion is also faster than v14 — roughly 3× on a large
alignment and 10–14× on the inputs that made v14 slowest — because the codon
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
| `-blockonly` | keep only columns marked `#` under the alignment |
| `-codontable N` | NCBI genetic code, default 1 |
| `-html` | HTML output |
| `-nostderr` | suppress messages |
| `-h` | usage |

Supported codon tables: 1–6, 9–16, 21–33 (every genetic code NCBI
defines), plus 34–37. NCBI assigns no ID past 33; 34–37 are the alternative
bacterial codes of Shulgina & Eddy (2021) under Wikipedia's numbering, which
a future NCBI assignment could collide with. Tables 7 and 8 were retired by
NCBI in 1995 and 17–20 have never been assigned, so neither is accepted.

Mark in-frame stop codons with `*` or `_`, and frame shifts with a digit
giving the number of nucleotides the column consumes.

## Tests

```sh
python -m venv .venv && .venv/bin/pip install -e '.[dev]'
.venv/bin/python -m pytest
```

The suite compares the port against output captured from the original Perl.
`tests/golden/` holds what `pal2nal.v14.pl` produces for every case in
`tests/cases.tsv`; `tests/expected/` holds the intended output for the cases
where v15 deliberately differs — 47 of the 145 cases, each registered in
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
| `CHANGELOG.md` | what changed in v15 and why |
| `PORTING.md` | the porting decisions and the evidence behind them |

The web front end, the full Perl and CGI history of the original EMBL
service, and the notes on how that history was reconstructed live in a
separate archive repository and are deliberately not part of this package.

## Licence

MIT — see [LICENSE](LICENSE).
