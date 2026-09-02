# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PAL2NAL converts a protein sequence alignment plus the corresponding DNA
sequences into a codon alignment. This repository is **v15**, a Python port
of `pal2nal.pl` v14 (Perl, 2011, Mikita Suyama). A copy of the original is
kept at `tests/reference/pal2nal.v14.pl` and is the reference implementation.

This is the command-line tool only. The Flask web front end and the full
reconstructed history of the original EMBL Perl/CGI service live in a
separate archive repository.

## Commands

```sh
python -m venv .venv && .venv/bin/pip install -e '.[dev]'

.venv/bin/python -m pytest                      # whole suite
.venv/bin/python -m pytest tests/test_golden.py # port vs. the Perl only
.venv/bin/python -m pytest -k dhfr_output_paml  # one case, by name

cd tests && sh generate_golden.sh               # recapture the Perl's output
cd tests && sh generate_golden.sh /path/to/pal2nal.v12.2.pl   # another version
```

Requires `perl` only for `generate_golden.sh`. There is no configured linter
or formatter; `pyproject.toml` sets ruff's line length to 100.

## The golden corpus is the specification

This is the thing to understand before changing anything.

`tests/golden/` holds byte-exact stdout, stderr and exit status captured
from `pal2nal.v14.pl` for every case in `tests/cases.tsv`. The port must
reproduce it exactly — wording, spacing, column widths, which stream a
message goes to, exit status. **Those files are evidence. Never edit one to
make a test pass, and never regenerate them to paper over a difference.**

Where v15 is *meant* to differ, the case is listed in
`tests/divergences.tsv` with a prose reason, and the intended output lives
in `tests/expected/`. `test_golden.py` then compares against `expected/`
instead of `golden/`, so the Perl's behaviour stays on record and the change
is a reviewable diff rather than a silently edited expectation. A guard test
fails if a divergence lacks a reason or an `expected/` file.

To add a divergence: make the change, append a line to
`tests/divergences.tsv`, run `.venv/bin/python tests/record_expected.py`,
then **read the recorded output and confirm it is correct**. That script
records whatever the port currently does; it proves nothing on its own.

Every behavioural difference from v14 must also appear in `CHANGELOG.md`.

## Architecture

The pipeline mirrors the Perl's control flow, one module per stage:

```
options.py  parse argv          → inputs.py  read + normalise
            → validate.py match the input against the expected alphabet
            → convert.py  match peptides to codons
            → output.py   reinsert gaps, filter, format
            → cli.py      orchestrate, route messages to stdout/stderr
```

- `codontables.py` — hand-transcribed from the Perl's `%p2c` hashes (there
  is no generator script), plus tables 24-37 added in v15: amino acid →
  *regex fragment* matching the codons that could encode it. Alternations
  spell out IUPAC ambiguity (`U|T`, `A|G|R`) and `.` is a degenerate third
  position, so ambiguity codes in the input still match. Two keys are not
  amino acids: `B` is the table's initiation codons, `X` is any codon.
- `convert.py` — the core. Concatenate every residue's fragment into one
  pattern and search the DNA unanchored and case-insensitively. If it
  matches, every codon is correct by construction and nothing is checked.
  If not, cut the peptide into ten-residue anchors, relax any anchor that
  cannot be found on its own to wildcards, match the mixed pattern, then
  verify each codon individually so mismatches can be reported.
- `output.py` — `pn2codon` returns codons with no gaps; gaps are put back by
  walking the peptide alignment column by column. Column width is set by the
  widest entry, so a frame shift consuming four nucleotides in one sequence
  widens that column for all of them and the rows stay in register.
- `validate.py` — the gates `inputs.py` calls before and after parsing.
  Non-ASCII input (undecodable bytes included, via `surrogateescape`) and
  control characters other than tab/CR/LF are refused with a count and a
  position and **never an echo of the input**; DNA outside the IUPAC
  alphabet is refused by name. Peptide residues stay lenient — unknown
  characters are taken as X and reported per position, as v14 did — but
  no output filter may hide that report. Error behaviour is explicitly not
  held to v14's; see "Agreed for v15" in `PORTING.md`.

Note that `pn2codon` is quadratic in sequence length when the peptides and
the DNA do not correspond. The command-line tool applies no time limit; any
caller exposing it to untrusted input has to bound it itself.

## Behaviours that look like bugs but are deliberate

Reproduced from v14 on purpose, because the goldens pin them:

- Options are parsed in one left-to-right pass. `-output` and `-codontable`
  consume the *next* token whatever it is.
- `Options.alnfile` is assigned under a truthiness test, matching Perl's
  `!$alnfile`, so a file literally named `0` behaves oddly.
- The CLUSTAL conservation/marker line is located by **position**, never by
  content: whatever line follows a sequence line is sliced at the same
  columns. Changing this risks altering `-blockonly` on real alignments.
- `perl_num()` in `options.py` mimics Perl's loose numeric coercion, so
  `-codontable 1abc` is accepted as 1.
- `perl_split_lines()` drops trailing empty fields, as Perl's `split` does.

Twelve genuine defects *were* fixed; `CHANGELOG.md` lists them with reasons
and `PORTING.md` records the decisions. Notably `-codontable 10` (Euplotid
nuclear) was never implemented in v14 and silently returned universal-code
results.

`tests/test_codontables.py` decodes every table back into a codon-to-residue
mapping and diffs it against NCBI's `gc.prt`, vendored as
`tests/data/ncbi_gc.prt`. That check is what found the table 10 defect, and
it is the first thing to run after touching a table. Tables 34-37 are **not**
NCBI assignments — the numbering is Wikipedia's; see `PROVENANCE` in
`codontables.py` before adding more.

## Gotchas

- `.gitattributes` sets `* -text`. `tests/data/dhfr_crlf.aln` and
  `dhfr_cr.aln` carry Windows and classic-Mac line endings deliberately;
  checkout normalisation would flatten them and silently void those cases.
  This has already happened once.
- Goldens embed the input paths from the command line, so `cases.tsv` uses
  relative paths and both the runner and `test_golden.py` must execute from
  `tests/`.
- dS/dN via PAML `codeml` and the `bl2seq` diagnostic were dropped; don't
  reintroduce them without asking.
