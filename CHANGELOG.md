# Changelog

## v15.0 (2 September 2026)

First release since v14 (2 December 2011), and the first in Python. The
conversion algorithm is unchanged: v15 reproduces v14's output byte for
byte on every input except where a defect is fixed below. Each fix has a
regression test, and v14's own output is kept alongside the new expectation
so every change is visible as a diff (`tests/golden/` vs `tests/expected/`).

### Added

* **Fourteen more genetic codes.** v14 covered every NCBI table that existed
  when it was written in 2011 and rejected anything else, so the ten codes
  NCBI has published since (24-33) were unavailable. All ten are now
  supported, derived from NCBI's `gc.prt` (vendored as
  `tests/data/ncbi_gc.prt`). In tables 27, 28 and 31 every stop codon also
  codes for an amino acid; those codons match both, as NCBI's data says.
* **Four non-NCBI codes, 34-37.** The alternative bacterial codes found by
  Shulgina & Eddy (2021, eLife 10:e71402), numbered 34-37 by Wikipedia's
  "List of genetic codes". NCBI has neither adopted them nor assigned any id
  past 33, so **these numbers are Wikipedia's and could collide with a future
  NCBI assignment**. No initiation codons are published for them; the
  bacterial set (table 11) is assumed. `PROVENANCE` and `ASSUMED_STARTS` in
  `pal2nal/codontables.py` record both caveats.
* **`tests/test_codontables.py`**, which decodes every table back into a
  codon-to-residue mapping and checks it against its published source on all
  64 codons. Running it against v14's own hashes reports exactly one wrong
  codon in the seventeen tables it shipped -- TGA in table 10, below -- which
  is how that defect was identified; the test keeps that result on record.
  It also pins the initiation codons, so a future gap cannot pass unnoticed.

### Added: input validation

v14 read whatever bytes it was handed. Nothing checked the alphabet: a
character that was not a residue became "X" with a warning, a stray letter
in the DNA was kept and quietly broke the codon match, and anything that
was not a letter at all was deleted in silence. When the match then failed,
the error dumped both sequences back at the caller -- under `-html`, into
the page it returned. v15 gates both input files first, in
`pal2nal/validate.py`.

* **Non-ASCII input is refused**, and the message says only how many
  offending characters there are and where the first one is. The input is
  never echoed. This catches the gap typed as an en dash, the UTF-8 byte
  order mark a Windows editor prepends -- v14 read it as part of the first
  `>` line, lost that sequence and blamed the sequence count -- and a
  latin-1 file, which v14 converted silently and wrote the undecodable byte
  straight into an output ID. Undecodable bytes arrive as surrogates and
  are refused with the rest, so a binary upload is rejected rather than
  pattern-matched.
* **ASCII control characters other than tab, CR and LF are refused**, named
  by code point rather than written out. v14 echoed an embedded NUL back to
  the terminal.
* **DNA characters outside the IUPAC alphabet are refused by name and
  position.** Protein pasted into the DNA field is the commonest way to
  reach v14's generic "inconsistency between the following pep and nuc
  seqs"; that message is now reached only by input that really is DNA. The
  report is capped at five sequences, so a wholly corrupt file produces a
  report and not a second data dump. Only letters are checked, because
  digits, dashes and whitespace are stripped by design -- numbered and
  gapped FASTA variants depend on it.
* **`-nomismatch` and `-blockonly` no longer hide a character that is not a
  residue.** Both select which *codons* are reported; in v14 they also
  silenced "unknown AA type", so a column could be dropped with no
  explanation. The alignment they produce is unchanged. `-nostderr`, which
  asks for silence outright, still suppresses everything.
* **An error under `-html` now closes the `<pre>` it opened**, so a
  rejected input cannot leave the page's markup unbalanced.

Peptide residues stay lenient: an unknown character is still taken as X and
reported per alignment position, as v14 did, because B, J, O and Z are real
ambiguity codes and the conversion is still useful. An unknown nucleotide
has no such fallback and is fatal.

### Fixed: wrong results

* **Codon table 10 (Euplotid nuclear) was never implemented.** TGA was left
  in the stop pattern and missing from Cys, so `-codontable 10` silently
  produced universal-code results. TGA now codes for Cys as NCBI defines
  it. Found by decoding all seventeen of v14's tables back into codon
  assignments and diffing them against the official NCBI translation
  strings: table 10's TGA is the only disagreement in 1088 codons, and the
  other sixteen tables are unchanged. `tests/test_codontables.py` reruns
  that comparison.
* **Three tables rejected legitimate alternative initiation codons.** v14's
  start-codon lists predate later NCBI revisions: the yeast mitochondrial
  code (3) did not accept GTG, the ascidian mitochondrial code (13) did not
  accept ATA — NCBI added both after v14 was written — and the bacterial
  code (11) omitted ATT, ATC and ATA, with no NCBI changelog entry to explain
  it. A gene with any of those starts drew a spurious "pepAlnPos 1: M does
  not correspond to ..." warning, which matters for bacteria, where GTG and
  TTG starts are common. All three now accept exactly what NCBI lists. Only
  the warning changes; the codon alignment itself was always correct.
* **Alignments of ten or fewer residues could not use the fallback matcher.**
  The step that merges a short trailing anchor iterated an empty range when
  there was only one anchor, leaving the search pattern undefined; any short
  alignment needing that path failed with "inconsistency between the
  following pep and nuc seqs". A five-residue alignment with one mismatched
  codon now converts, as the same mismatch at twelve residues always did.
* **Selenocysteine aborted the run.** `U` was accepted as a residue but no
  codon table defined it, so it contributed an empty pattern that
  desynchronised everything after it. `U` now matches TGA.
* **A protein alignment in FASTA format beginning with a blank line was
  parsed as CLUSTAL.** The format was detected from the file read as a
  single string, so the leading newline defeated the test and the alignment
  was split on whitespace into as many "sequences" as it had words.
  Detection now starts at the first non-blank, non-comment line.
* **`.` was documented as a gap but only honoured as one while matching**,
  not while writing the alignment, so a `.` column produced no gap where a
  `-` column did. Both are now treated alike throughout.
* **The fallback matcher used a fixed set of initiation codons** rather than
  the selected table's own, and that fixed set matches none of the seventeen
  tables. Each table's initiation codons are now used.
* **`-output codon` could misalign the amino-acid row.** The width of a
  widened frame-shift column was taken from a leftover loop variable
  instead of the column maximum, so when the frame-shift digit belonged to
  any sequence but the last, the amino-acid row came out a column short and
  drifted against the codon row from that point on.
* **A repeated FASTA header merged two records but counted them twice**,
  emitting the merged sequence for both. Duplicate IDs are now rejected.
  CLUSTAL input already deduplicated correctly.
* **ID-based pairing could be chosen when the IDs did not correspond.** The
  choice between ID-based and positional pairing compared list lengths, not
  the ID sets, so duplicates could satisfy it and a protein would be paired
  with a missing DNA sequence. A genuine one-to-one match is now required.

### Changed: command-line behaviour

* **Errors exit with status 1.** Every validation failure previously exited
  0, so `pal2nal ... || handle_error` never fired and pipelines silently
  continued on a failed conversion.
* **`-h` prints usage and stops.** It previously printed usage and then ran
  the conversion anyway.
* **A missing option value is an error.** A trailing `-output` or
  `-codontable` was ignored and the default used, and `-output -html`
  consumed `-html` as the format name.

### Removed

* **dS/dN calculation.** It depended on a 32-bit Linux `codeml` binary from
  around 2004 and parsed a PAML output format that current releases no
  longer produce. The tree file it generated was empty in any case: the
  variable holding the sequence names was never populated, so every run
  wrote `(, );`.
* **The `bl2seq` diagnostic** printed on a protein/DNA inconsistency. The
  tool is retired; the inconsistency message itself is unchanged.
* **The mail-a-copy route** (`sendcopy.cgi`, `sendcopy.pl`), which shelled
  out to `mutt` and was pinned to the v8 tarball.
* **`clean_tmp.pl`**, the cron job that swept the CGI's scratch files.

### Notes

* Not fixed, deliberately: the CLUSTAL conservation line is still located by
  position rather than content, because changing that risks altering
  `-blockonly` results on real alignments.
* Error and diagnostic output is **not** held to v14's wording or exit
  status. What v14 did with each rejected input is still captured in
  `tests/golden/` and what v15 does in `tests/expected/`, so every change
  remains a reviewable diff.
* `tests/reference/pal2nal.v14.pl` keeps the original Perl; `PORTING.md`
  records the porting decisions. The full Perl and CGI history of the EMBL
  service, and the notes reconstructing it, are in a separate archive
  repository.

## v14 and earlier

See the changelog block at the top of `tests/reference/pal2nal.v14.pl`. The
earlier releases, replayed commit by commit, are in the separate archive
repository; they were kept under version-stamped file names rather than as
revisions, so reading their diffs needs `git log -C --find-copies-harder`.
