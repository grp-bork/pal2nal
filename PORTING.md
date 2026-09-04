# Porting pal2nal from Perl to Python

v15 is a Python port of `pal2nal.pl` v14 (2 December 2011), the last version
deployed on the EMBL web server. The Perl is vendored at
`tests/reference/pal2nal.v14.pl` and was the reference throughout: v15
reproduces its output byte for byte except where this file records a
deliberate divergence.

This document records the decisions the port made and the evidence behind
them. `CHANGELOG.md` lists the resulting user-visible changes;
`tests/README.md` describes the corpus.

## How the port was verified

* `tests/golden/` holds the Perl's own output for every case in
  `tests/cases.tsv` (`wc -l < tests/cases.tsv`), captured by
  `tests/generate_golden.sh`. Those files are evidence and are never edited
  to make a test pass.
* Where the port deliberately differs, the case is listed in
  `tests/divergences.tsv` with a reason and the intended output is stored
  in `tests/expected/` — 75 of the 171 cases as of this writing; both
  counts grow as options are added, so `wc -l < tests/cases.tsv` and the
  non-comment, non-blank lines of `tests/divergences.tsv` are the numbers
  to trust over this prose. The Perl's behaviour stays on record in
  `golden/`, so each difference is a reviewable file rather than a silent
  edit.
* Everything else about a *successful conversion* — wording, spacing,
  column widths, which stream a message goes to — was treated as specified
  behaviour.
* Error behaviour was not. Diagnostics, wording and exit status changed
  where v14's were wrong or unhelpful; every such change is recorded too,
  so what v14 did stays visible.

## Genetic codes

### Table 10 (Euplotid nuclear) was never implemented

The Euplotid nuclear code deviates from the universal code in exactly one
codon: **TGA codes for Cys** instead of terminating translation.
`pal2nal.v14.pl` does not implement this. In its table 10 (line 1450):

    "C" => "((U|T)G(U|T|C|Y))"                    # TGT, TGC only
    "_" => "(((U|T)A(A|G|R))|((T|U)GA))"          # TGA still a stop
    "*" => "(((U|T)A(A|G|R))|((T|U)GA))"

TGA is absent from the Cys pattern and still present in both stop
patterns, so `-codontable 10` produced results identical to
`-codontable 1`. The test corpus confirms it in both directions:
`codontable_10_own_probe` reports a mismatch on every TGA column, and
`codontable_t1_under_10` converts silently where it should not.

The port implements the table as NCBI defines it. Anyone who used
`-codontable 10` with v14 or earlier received universal-code results, so
output for Euplotid data changes — correctly.

The other 16 tables v14 shipped (1, 2, 3, 4, 5, 6, 9, 11, 12, 13, 14, 15,
16, 21, 22, 23) agree with the official NCBI `ncbieaa` strings on all 64
codons and were ported unchanged.

How this was established, and how it stays established: every table is
decoded back from its regexes into a plain codon-to-residue mapping and
diffed against NCBI's `gc.prt`. Over v14's seventeen tables that is 1088
codons, and table 10's TGA is the single disagreement.
`tests/test_codontables.py` runs that comparison against both the port and
the Perl, so the finding is reproducible rather than a claim in prose.

### Fourteen genetic codes were added

v14 shipped every NCBI table that existed in 2011. Tables 24–33 have been
published since and were added from `gc.prt`; `-codontable 24`..`33` were
rejected as invalid by v14, so each has a divergence entry. Tables 34–37
are the four alternative bacterial codes of Shulgina & Eddy (2021), which
NCBI has *not* adopted: the numbering is Wikipedia's and is not
authoritative. See `PROVENANCE` in `pal2nal/codontables.py`.

Two details are worth knowing. NCBI marks a context-dependent stop in the
`sncbieaa` row rather than `ncbieaa`, so in tables 27, 28 and 31 a codon can
be both a residue and a stop; those codons match both patterns. And no
initiation codons are published for 34–37, so the bacterial set from table
11 is assumed and listed in `ASSUMED_STARTS`.

### v14's initiation codons were widened to match NCBI

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

## Defects fixed

Sixteen genuine defects were fixed. The numbers below are the ones
`tests/divergences.tsv` cites; `CHANGELOG.md` describes the same sixteen for
users. In each case v14's behaviour stays recorded in `tests/golden/` and
the new behaviour in `tests/expected/`.

### Silent wrong answers

1. **Alignment format detection.** v14 sniffs the format from the file read
   as a single string (`undef $/` at line 242 turns the `while (<ALNFILE>)`
   at line 366 into one whole-file iteration), so `^CLUSTAL` and `^>` only
   match at byte 0. A FASTA alignment beginning with a blank line was
   parsed as CLUSTAL and fell apart — reproduced: two sequences became
   "aa: 9". v15 detects from the first non-blank, non-comment line.
2. **Duplicate FASTA peptide IDs.** v14 appended a repeated `>id` to the ID
   list twice while merging both records' residues under one key, emitting
   the merged sequence twice and inflating the sequence count. CLUSTAL
   input deduplicated properly. v15 rejects duplicate IDs with an error.
3. **ID-based vs positional pairing.** v14 chose between them by comparing
   list *lengths* (`$#commonids == $#aaid`, line 488), which duplicates can
   satisfy without the ID sets matching, silently pairing a peptide with a
   missing DNA sequence. v15 requires a genuine one-to-one match.

### Command-line contract

4. **Exit status.** Every validation failure in v14 exited 0, so
   `pal2nal ... || handle_error` never fired; only a failed file open
   exited non-zero. v15 exits 1 on error.
5. **`-h`.** v14 printed usage and then carried on converting when two or
   more arguments were present. v15 prints usage and stops.
6. **Missing option values.** `-output` or `-codontable` as the final token
   was silently ignored and the default used, and `-output -html` consumed
   `-html` as the format. v15 reports an error.
15. **A non-existent option in an error message.** Rejecting `-output
   codon` alongside a filter, v14 reported `"-outform codon" is not valid
   with -blockonly, -nogap, -nomismatch`. The tool has never had an
   `-outform` option, so following the message produced "invalid output
   format" from the next run. v15 names `-output`.

### Matching and output

7. **The short-peptide fallback.** The step that merges a short trailing
   anchor iterated an empty range when there was only one anchor, leaving
   the search pattern undefined, so any alignment of ten or fewer residues
   that needed the fallback aborted with "inconsistency between the
   following pep and nuc seqs". v15 converts it and reports the mismatched
   codon, as the same mismatch at twelve residues always did.
8. **Selenocysteine.** `U` was accepted as a residue but no codon table
   defined it, so it contributed an empty pattern that desynchronised
   everything after it and aborted the run. v15 matches `U` to TGA.
9. **Codon-row drift under `-output codon`.** The width of a widened
   frame-shift column was taken from a leftover loop variable instead of
   the column maximum, so when the frame-shift digit belonged to any
   sequence but the last, the amino-acid row came out a column short and
   drifted against the codon row from that point on. v15 uses the column
   maximum.
10. **Initiation codons in the fallback matcher.** The fallback path used a
    fixed set of initiation codons rather than the selected table's own,
    and that fixed set matches none of the seventeen tables v14 shipped.
    v15 uses each table's own.
11. **`.` as a gap.** It was documented as a gap and honoured as one while
    matching, but not while writing the alignment, so a `.` column produced
    no gap where a `-` column did. v15 treats both alike throughout.
12. **HTML escaping.** v14 wrote IDs and sequences into `-html` output raw,
    so an ID containing `<` or `&` produced broken markup. v15 escapes
    them, and an error under `-html` closes the `<pre>` it opened.
13. **`-nogap` stop codons.** The filter tested for a stop with the
    universal code's TAA/TAG/TGA written out literally (line 857), ignoring
    `-codontable`, so on every other code it erred in both directions at
    once: under table 6 it deleted the TAA and TAG columns that spell Gln,
    and under table 2 it deleted the TGA column that spells Trp while
    keeping the AGA and AGG stops the option exists to remove. `-nogap` is
    what prepares an alignment for dS/dN, so this both dropped good codons
    and let a real in-frame stop through. v15 uses the selected table's own
    `*` pattern, the same one the matcher uses.
16. **Ragged alignments.** The alignment length was read off whichever
    row happened to be first, so a file whose rows differed in width was
    mangled one of two ways depending only on the order the sequences
    appeared in. With the long row first the short ones simply ran out and
    printed a row a column short; with the short row first the tail of
    every longer row was dropped, taking real codons with it. Both exited
    0 with no message. v15 refuses the file, names the rows and their
    lengths, and does not guess which width was meant.
14. **The Gblocks parser.** `$getaln` was reset to 0 at the top of every
    line of a Gblocks-format file, before the branch that consumes the
    data (line 433), so every sequence line was discarded and the run died
    with "number of input seqs differ (aa: 0; nuc: 2)". The format could
    never have worked. v15 keeps the flag across lines; the `#` mask on the
    Gblocks line then also reaches `-blockonly`.

## Input validation was hardened

Error behaviour is explicitly **not** held to v14's. The goldens still
record what v14 did with each of these inputs, which is the evidence for
why the behaviour changed.

* **Non-ASCII and control characters are refused** by `validate.py`,
  before parsing, and the message reports a count and a position rather
  than the offending text. v14 had no encoding contract at all: it matched
  the UTF-8 bytes of an en dash as residues and echoed the sequence back
  on failure (into the page under `-html`); it read a UTF-8 BOM as part of
  the first `>` line, lost that sequence and reported a count mismatch;
  and it converted a latin-1 file silently, writing the undecodable byte
  into an output ID. Files are read with `errors="surrogateescape"`, so
  undecodable bytes are non-ASCII here too and a binary upload is rejected
  rather than pattern-matched.
* **DNA outside the IUPAC alphabet is refused by name.** A stray letter
  used to survive the `[^a-zA-Z]` strip, break the match, and produce the
  generic "inconsistency between the following pep and nuc seqs" dump —
  the usual way of learning that protein had been pasted into the DNA
  field. Only letters are checked: digits, dashes and whitespace are
  stripped deliberately, because numbered and gapped FASTA variants
  depend on it. Peptide residues stay lenient by contrast, still taken as
  X with a per-position report, because B, J, O and Z are real ambiguity
  codes and the conversion remains useful.
* **`-nomismatch` and `-blockonly` no longer suppress that report.** They
  choose which codons are reported; neither can make a character that is
  not a residue acceptable. The alignment they emit is byte-identical to
  v14's — only the report is added. `-nostderr` still silences
  everything. `convert.py` returns the two kinds of message in separate
  fields so the distinction is structural rather than a string test.

## `-partial`: recovering a codon alignment from a broken pep/CDS pair

New in v16, and not a port of anything v14 did, so it gets its own section.

v14 and v15 abort the entire run on the first peptide/CDS pair that does not
correspond, with "inconsistency between the following pep and nuc seqs" —
discarding every other sequence in the alignment along with the offending
one. Real Ensembl and NCBI data hits this routinely, badly enough that a
third-party fork, `dukecomeback/pal4nal.pl`, exists solely to work around
it, by shelling out to `genewise`. `-partial` solves the same problem
in-process: no external aligner, no dependency to install, and it reports
what it could not place rather than guessing at it.

`pn2codon` returns `NO_MATCH` for exactly four shapes of input: an intron
in the DNA, an indel of any size between the peptide and its CDS, and a CDS
truncated at either end. (A point substitution already worked — it returns
`MISMATCH`, not `NO_MATCH`, because v14's relaxed-wildcard fallback finds
it.) All four share one property: the DNA between two of the peptide's
ten-residue anchors is not the width the anchors assume. `_chain()` in
`pal2nal/convert.py` places each anchor independently, at the leftmost DNA
offset at or after the end of the previously placed anchor, instead of
concatenating every anchor into one fixed-width pattern — so the gap
between two anchors can be any length, which is what all four pathologies
need. It uses `Offsets.matches()`, split out of `Offsets.search()` in
`pal2nal/codonmatch.py`, which returns the whole surviving offset bitmask
rather than collapsing it to the lowest set bit.

**The chain is gated on `NO_MATCH`.** It runs only after both of v14's
matching paths — the whole-pattern search and the relaxed-anchor fallback —
have already failed. This was checked empirically, not assumed: chaining
unconditionally is worse than v14's fallback on cases that already convert.
The decisive case is `mismatch_warned`, sequence `m2`: peptide `MAKQLRT`
against `ATGGCTAAAGGGCTGCGTACT`, one substitution. v14's fallback recovers
all seven codons and reports the one mismatch. Chaining recovers nothing,
because a single seven-residue anchor carrying a substitution cannot be
found on its own and has no neighbouring anchor to bracket it. Because of
the gate, `-partial` provably changes nothing about any input that
converts today — `tests/test_partial.py` checks this directly, by running
every non-`partial_*` case in the corpus twice, with and without the flag,
and asserting the codon alignment and exit status come out identical.

`_fill_runs()` then recovers an unplaced anchor when it is bracketed by two
placed anchors and the DNA span between them is exactly the run's width —
this is what keeps a point substitution like `m2` from being worse under
chaining than under the old fallback. The fill applies only to interior
runs: a run at either end of the sequence is pinned on one side only, and
filling it would read into the 5' UTR or the 3' UTR as if it were coding
sequence, fabricating a codon from whatever happens to flank the peptide
rather than reporting that nothing is there.

Anchors are chained **greedily**, left to right. A tandem repeat can attract
an anchor to the wrong copy of itself; this is a known limitation, not a
bug to be fixed later, and the residue counts in the `-partial` report
(`N/M residues placed`) are what would reveal it on real data.

The gap filler for an anchor that cannot be placed is `"-" * width`, where
`width` is computed the same way `_build_steps` computes it, never three
dashes per residue. A frame-shift numeral consumes `int(aa)` nucleotides,
so for the peptide `EK4QKNDTY` the true width is 28, not `3 * 9 = 27`.
Under-spending by even one nucleotide there shifts every downstream codon
by one base — `tests/data/drift_pep.fasta` pins that same defect one level
down, independent of `-partial`, and `partial_frameshift` is the case that
exercises it here.

`_build_steps` gained a `seen_letter` keyword so that `_chain` can thread
it across anchors as they are placed one at a time. Without it, every
anchor beginning with M would match the table's *initiation* codons — the
`B` pattern — rather than the ordinary Met codon, because `_build_steps`
used to reset its own notion of "have we seen the first residue yet" on
every call. Under table 1, `B` is `((U|T|C|Y|A)(U|T)G)`, which also accepts
CTG and TTG; under table 11 it additionally accepts ATT, ATC, ATA and GTG.
An untethered `_chain` would let a mid-sequence Met match any of those,
silently, as if it had been verified.

**An unrelated CDS is deliberately left alone.** It reaches v14's
relaxed-wildcard fallback, "succeeds" by matching all wildcards at offset
zero, and yields confident garbage — the codon alignment is emitted as if
it were correct, with nothing to flag it. That is v14 behaviour pinned by
the goldens, and it returns `MISMATCH` rather than `NO_MATCH`, so
`-partial`'s gate never lets the chain engage. Separating this case from a
legitimate one like `m2` would need a "fraction of codons that verify"
threshold, and the corpus discipline this project holds itself to should
not accept a threshold like that without a principled definition of where
it sits. Instead, under `-partial` the report gains a
`MISMATCH: <id> 1/60 codons verified` line for exactly this case, which is
what makes the result visibly worthless without pretending to distinguish
it algorithmically from a real one.

Finally: `-partial` relaxes only the peptide/CDS correspondence check.
Ragged alignments, duplicate IDs, bad alphabets and every other validation
gate still abort the run exactly as they did before. "Does not abort" is
not "accepts anything". Exit status stays 0 when `-partial` recovers
something, because not aborting is the point of the flag.

## Features not ported

* **`bl2seq`.** On a protein/DNA inconsistency v12 and later shelled out to
  the retired NCBI `bl2seq` to print an alignment as a diagnostic. The
  plain inconsistency message is kept; the alignment is not.
* **`codeml` / the dS-dN calculation.** The bundled `codeml` is a 32-bit
  Linux ELF binary from around 2004 and the CGI parsed a 2006-era output
  format. The conversion tool is ported; the Ks/Ka option is not.

The web front end is not part of this package either, and neither are the
pieces that only ever served it: `sendcopy.pl` / `sendcopy.cgi` /
`mail.txt`, which mailed the tarball via `mutt` from a hardcoded path and
were pinned to v8; `clean_tmp.pl`, the cron reaper for CGI scratch files;
the `rsh batman`, `wc`, `which` and `ls` shell-outs, replaced by in-process
equivalents; and the Matomo tracker in `index.cgi`, which reported to
`tr-denbi.embl.de`. They live in the separate archive repository along with
the rest of the EMBL service's history.

## Behaviour confirmed against the goldens

Everything here is v14 behaviour the port keeps, each one pinned by a case:

* Line endings: LF, CRLF and bare-CR input all produce identical output
  (v11/v13 behaviour). Covered by `aln_crlf_windows` and
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
* The CLUSTAL conservation line is located by position rather than
  content. This is a defect, deliberately not fixed: changing it risks
  altering `-blockonly` results on real alignments.
