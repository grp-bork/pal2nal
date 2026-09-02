# Porting pal2nal from Perl to Python

The reference implementation is `tests/reference/pal2nal.v14.pl` (2 December
2011), the last version deployed on the EMBL web server. The port must
reproduce it byte for byte, except where a divergence is recorded below.

## Ground rules

* `tests/golden/` holds the Perl's own output for all 127 cases in
  `tests/cases.tsv`, captured by `tests/generate_golden.sh`. Those files
  are evidence and are never edited to make a test pass.
* Where the port is meant to differ, the case is listed in
  `tests/divergences.tsv` with a reason and the intended output is stored
  in `tests/expected/`. The Perl's behaviour stays on record in
  `golden/`, and the difference is a reviewable file rather than a silent
  edit.
* Everything else about a *successful conversion* — wording, spacing,
  column widths, which stream a message goes to — is treated as specified
  behaviour.
* Error behaviour is not. Diagnostics, wording and exit status may change
  where v14's were wrong or unhelpful; the divergence is still recorded so
  what v14 did stays visible.

## Deliberate divergences

### `-codontable 10` (Euplotid nuclear) is fixed

The Euplotid nuclear code deviates from the universal code in exactly one
codon: **TGA codes for Cys** instead of terminating translation.
`pal2nal.v14.pl` does not implement this. In its table 10 (line 1450):

    "C" => "((U|T)G(U|T|C|Y))"                    # TGT, TGC only
    "_" => "(((U|T)A(A|G|R))|((T|U)GA))"          # TGA still a stop
    "*" => "(((U|T)A(A|G|R))|((T|U)GA))"

TGA is absent from the Cys pattern and still present in both stop
patterns, so `-codontable 10` produces results identical to
`-codontable 1`. Confirmed in both directions by the test corpus:
`codontable_10_own_probe` reports a mismatch on every TGA column, and
`codontable_t1_under_10` converts silently when it should not.

The port implements the table as NCBI defines it. Anyone who used
`-codontable 10` with v14 or earlier received universal-code results, so
output for Euplotid data will change — correctly.

The other 16 tables (1, 2, 3, 4, 5, 6, 9, 11, 12, 13, 14, 15, 16, 21, 22,
23) agree with the official NCBI `ncbieaa` strings on all 64 codons and
are ported unchanged.

How this was established, and how it stays established: every table is
decoded back from its regexes into a plain codon-to-residue mapping and
diffed against NCBI's `gc.prt`. Over v14's seventeen tables that is 1088
codons, and table 10's TGA is the single disagreement.
`tests/test_codontables.py` runs that comparison against both the port and
the Perl, so the finding is reproducible rather than a claim in prose.

### Genetic codes added in v15

v14 shipped every NCBI table that existed in 2011. Tables 24-33 have been
published since and are added from `gc.prt`; `-codontable 24`..`33` were
rejected as invalid by v14, so each has a divergence entry. Tables 34-37
are the four alternative bacterial codes of Shulgina & Eddy (2021), which
NCBI has *not* adopted: the numbering is Wikipedia's and is not
authoritative. See `PROVENANCE` in `pal2nal/codontables.py`.

Two details are worth knowing. NCBI marks a context-dependent stop in the
`sncbieaa` row rather than `ncbieaa`, so in tables 27, 28 and 31 a codon can
be both a residue and a stop; those codons match both patterns. And no
initiation codons are published for 34-37, so the bacterial set from table
11 is assumed and listed in `ASSUMED_STARTS`.

### v14's initiation codons are widened to match NCBI

Three tables accepted fewer initiation codons than NCBI lists: table 3 was
missing GTG, table 13 ATA (both added to `gc.prt` after v14 was written) and
table 11 ATT/ATC/ATA (absent from v14 with no NCBI changelog entry to
explain it). All three now accept exactly NCBI's set.

The `B` pattern applies only to the first residue of a sequence, and only
when it is Met, so the effect is narrow: a gene beginning at one of those
codons drew a spurious `pepAlnPos 1: M does not correspond to ...` warning.
The codon alignment was unaffected either way — verified by diffing the
output before and after on a GTG-, ATT- and ATA-initiated sequence, where
only the warning disappears. No golden case exercises this path, so no
recorded expectation changed.

`test_start_codons_widened_exactly_where_recorded` decodes v14's own `B`
patterns and asserts the port adds exactly these codons and removes none;
`test_initiation_codons_match_the_published_table` holds every table to
NCBI's list.

## Dropped features

* **`bl2seq`.** On a protein/DNA inconsistency v12+ shelled out to the
  retired NCBI `bl2seq` to print an alignment as a diagnostic. The plain
  inconsistency message is kept; the alignment is not. No BLAST in the
  image.
* **`codeml` / the dS-dN calculation.** The bundled `codeml` is a 32-bit
  Linux ELF binary from ~2004 and the CGI parses a 2006-era output
  format. The conversion tool and the web front end are ported; the
  Ks/Ka option is not.
* **`sendcopy.pl` / `sendcopy.cgi` / `mail.txt`.** Mailed the tarball via
  `mutt` from a hardcoded path, pinned to v8.
* **`clean_tmp.pl`.** Cron reaper for CGI scratch files, unnecessary once
  temporary files are handled properly.
* **`rsh batman`, `wc`, `which`, `ls` shell-outs.** Replaced by in-process
  equivalents.
* **The Matomo tracker** in `index.cgi`, which reported to
  `tr-denbi.embl.de`.

## Things the port must get right

Recorded here as they are confirmed against the goldens:

* Line endings: LF, CRLF and bare-CR input must all produce identical
  output (v11/v13 behaviour). Covered by `aln_crlf_windows` and
  `aln_cr_classicmac`, which hash-match `dhfr_default_format`.
* Peptide input format is auto-detected; `-a` has not existed since v8.
* DNA may be upper or lower case (v10) and may be split across several
  files (v2).
* Sequences are paired by ID when the peptide and DNA IDs match, and by
  input order otherwise (v12).
* `-output codon` is rejected together with `-nogap`, `-blockonly` or
  `-nomismatch`.
* Warnings go to STDERR normally, to STDOUT under `-html`, and are
  suppressed under `-nostderr`.

## Agreed for v15

The port is released as **v15**, and every divergence below is a changelog
entry. In each case the Perl's behaviour stays recorded in `tests/golden/`
and the new behaviour in `tests/expected/`.

### Fixed: silent wrong answers

1. **Alignment format detection.** v14 sniffs the format from the file read
   as a single string (`undef $/` at line 242 turns the `while (<ALNFILE>)`
   at line 366 into one whole-file iteration), so `^CLUSTAL` and `^>` only
   match at byte 0. A FASTA alignment beginning with a blank line was
   parsed as CLUSTAL and fell apart -- reproduced: two sequences became
   "aa: 9". v15 detects from the first non-blank, non-comment line.
2. **Duplicate FASTA peptide IDs.** v14 appended a repeated `>id` to the id
   list twice while merging both records' residues under one key, emitting
   the merged sequence twice and inflating the sequence count. CLUSTAL
   input deduplicated properly. v15 rejects duplicate ids with an error.
3. **ID-based vs positional pairing.** v14 chose between them by comparing
   list *lengths* (`$#commonids == $#aaid`, line 488), which duplicates can
   satisfy without the id sets matching, silently pairing a peptide with a
   missing DNA sequence. v15 requires a genuine one-to-one match.

### Fixed: command-line contract

4. **Exit status.** Every validation failure in v14 exited 0, so
   `pal2nal ... || handle_error` never fired; only a failed file open
   exited non-zero. v15 exits 1 on error.
5. **`-h`.** v14 printed usage and then carried on converting when two or
   more arguments were present. v15 prints usage and stops.
6. **Missing option values.** `-output` or `-codontable` as the final token
   was silently ignored and the default used, and `-output -html` consumed
   `-html` as the format. v15 reports an error.

### Hardened: what the converter refuses

Error behaviour is explicitly **not** held to v14's. The goldens still
record what v14 did with each of these inputs, which is the evidence for
why the behaviour changed.

7. **Non-ASCII and control characters are refused** by `validate.py`,
   before parsing, and the message reports a count and a position rather
   than the offending text. v14 had no encoding contract at all: it matched
   the UTF-8 bytes of an en dash as residues and echoed the sequence back
   on failure (into the page under `-html`); it read a UTF-8 BOM as part of
   the first `>` line, lost that sequence and reported a count mismatch;
   and it converted a latin-1 file silently, writing the undecodable byte
   into an output ID. Files are read with `errors="surrogateescape"`, so
   undecodable bytes are non-ASCII here too and a binary upload is rejected
   rather than pattern-matched.
8. **DNA outside the IUPAC alphabet is refused by name.** A stray letter
   used to survive the `[^a-zA-Z]` strip, break the match, and produce the
   generic "inconsistency between the following pep and nuc seqs" dump --
   the usual way of learning that protein had been pasted into the DNA
   field. Only letters are checked: digits, dashes and whitespace are
   stripped deliberately, because numbered and gapped FASTA variants
   depend on it. Peptide residues stay lenient by contrast, still taken as
   X with a per-position report, because B, J, O and Z are real ambiguity
   codes and the conversion remains useful.
9. **`-nomismatch` and `-blockonly` no longer suppress that report.** They
   choose which codons are reported; neither can make a character that is
   not a residue acceptable. The alignment they emit is byte-identical to
   v14's -- only the report is added. `-nostderr` still silences
   everything. `convert.py` returns the two kinds of message in separate
   fields so the distinction is structural rather than a string test.
