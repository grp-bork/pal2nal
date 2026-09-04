# Test corpus

Golden outputs captured from the original implementation, vendored here
as `reference/pal2nal.v14.pl`, running under Perl 5.34. They are the
reference the Python port reproduces.

    sh generate_golden.sh                      # regenerate from v14
    sh generate_golden.sh /path/to/pal2nal.v12.2.pl   # or any other version

The goldens were originally captured with the interpreter at
`../legacy_code/pal2nal.v14.pl`, and `golden/err_missing_file.err` embeds
that path because Perl's `die` prints the script name and line number.
Regenerating from `reference/pal2nal.v14.pl` therefore rewrites that one
file. The case is a registered divergence, so `test_golden.py` compares it
against `expected/`, not against the golden.

`cases.tsv` lists each case as `name<TAB>arguments`. For every case the
runner records `golden/<name>.out` (stdout), `golden/<name>.err` (stderr,
omitted when empty) and `golden/<name>.status` (exit status). 171 cases
(`wc -l < cases.tsv`).

The script echoes its input file names into the output header, so the
paths in `cases.tsv` are relative and the runner must be started from
this directory.

## Where v15 differs

75 of the 171 cases are registered in `divergences.tsv`, each with a prose
reason (non-blank, non-comment lines). For those, `test_golden.py` compares
against `expected/` instead of `golden/`, so v14's behaviour stays on
record and the change is a reviewable diff. A guard test fails if a
divergence lacks a reason or an `expected/` file. `PORTING.md` explains the
decisions; `CHANGELOG.md` describes them for users.

**The files in `golden/` and `expected/` are evidence: never edit one to
make a test pass, and never regenerate them to paper over a difference.**

## What is covered

| Group | Cases |
|---|---|
| Output formats | `clustal`, `paml`, `fasta`, `codon`, and the default |
| Filters | `-nogap`, `-nomismatch`, `-blockonly` and their combinations |
| Web modes | `-html`, `-nostderr` and combinations |
| Input handling | CLUSTAL and FASTA peptide input, LF / CRLF / CR line endings, DNA split across several files, lower-case DNA, ID-matched and ID-matched-but-reordered DNA |
| Codon tables | all 17 tables v14 shipped, each probed with all 64 codons, plus cross-checks in both directions against the universal code; one probe for each of the 14 codes v15 adds |
| Biology | in-frame stop codons, all frame-shift notations from `howto.html`, UTRs and polyA tails, protein/DNA mismatch (warned and removed), clean baseline |
| Diagnostics | usage with no arguments, `-h`, invalid output format, invalid codon table, the `-output codon` conflicts, missing input file |
| Input validation | a gap typed as an en dash, a UTF-8 BOM, a latin-1 file, an embedded NUL, protein pasted into the DNA field, an unknown residue on its own and under `-nomismatch` and `-blockonly` |
| `-partial` | an intron, an indel, a 5' truncation, a 3' truncation, an unrelated CDS, a frame shift inside an unplaceable anchor, an all-gap row, and one unmatchable sequence alongside two good ones, combined with `-nogap`, `-nomismatch`, `-nostderr`, `-html` and `-output codon` |

## Codon table probes

`data/allcodons.nuc` holds all 64 codons in NCBI order. For each table N,
`data/allcodons_t<N>_pep.fasta` is that table's translation of those 64
codons, taken from the official NCBI `ncbieaa` strings rather than from
the Perl, so the probe is an independent reference.

Running table N's probe under `-codontable N` should convert silently.
Under v14, 16 of its 17 tables do; under v15, all 31 do.

### Table 10 in v14

Table 10 is the Euplotid nuclear code, whose single deviation from the
universal code is **TGA = Cys**. `pal2nal.v14.pl` does not implement it:

    "C" => "((U|T)G(U|T|C|Y))"                    # TGT, TGC — TGA missing
    "*" => "(((U|T)A(A|G|R))|((T|U)GA))"          # TGA still a stop

so `-codontable 10` behaves exactly like `-codontable 1`. Two cases show
it: `codontable_10_own_probe` warns on both TGA columns, and
`codontable_t1_under_10` is silent where it should have complained.

The goldens record that behaviour as it is. v15 implements the table as
NCBI defines it, so both cases are registered divergences; see
`PORTING.md`.

## `-partial` fixtures

`data/` carries one pair of files per pathology `-partial` recovers from —
`<name>_pep.fasta` and `<name>.nuc`, each holding one sequence — plus one
alignment exercising the flag's actual point, that one bad pair no longer
takes the rest of the run down with it.

| Fixture | What it is |
|---|---|
| `intron` | The CDS has a 60 nt intron spliced in; the true CDS is `trunc3`'s DNA (its first 45 nt) followed by `trunc5`'s |
| `indel` | A single nucleotide deleted from an otherwise correct CDS |
| `trunc5` | The CDS is missing its first 15 codons |
| `trunc3` | The CDS is missing its last 15 codons |
| `unrelated` | DNA with no relationship to the peptide at all, kept to show `-partial` does not touch it — see below |
| `onebad` / `onebad_good_only` | A three-sequence alignment where the middle sequence cannot be matched, and the same alignment with that sequence removed, so the two survivors can be diffed byte for byte between the two runs |
| `partial_frameshift` | A frame-shift numeral inside an anchor that cannot be placed, followed by one that can, pinning that the gap filler spends the numeral's own width and not three nucleotides per residue |

`unrelated` reaches v14's relaxed-wildcard fallback and "succeeds" there,
so it returns `MISMATCH` rather than `NO_MATCH` and `-partial`'s gate never
lets it engage — `partial_unrelated` in `cases.tsv` pins that the codon
output is unchanged and that the report instead grows a
`MISMATCH: <id> 1/60 codons verified` line, which is what makes the result
visibly worthless.

`cases.tsv` pairs each fixture with a `partial_*` case running it under
`-partial` and, for the ones that would otherwise be silently mistaken for
success, a matching case without the flag (`partial_intron_off`,
`partial_multi_off`) so the abort v14 and v15 already produce stays on
record too. `partial_allgap` reuses the existing `allgap` fixture: an
all-gap row has nothing to place, so under `-partial` it converts straight
to the all-gap codon row it always implied instead of aborting.

## `tests/test_partial.py`

A free-standing test module for properties the golden corpus cannot express
by itself: that `-partial` is inert on everything that already converts,
that the codon string always keeps the alignment's register, and that the
gate and the two placement hazards (greedy chaining, initiation codons
threaded across anchors) behave as designed rather than by accident.

Its most valuable test, `test_partial_never_changes_an_alignment_that_
already_converted`, runs every non-`partial_*` case in `cases.tsv` twice,
once as recorded and once with `-partial` appended, and asserts the codon
alignment and exit status are identical between the two runs — the report
may gain a coverage line, but nothing already there may change. That is
the gate in `PORTING.md` checked directly against the whole corpus, rather
than argued for from one hand-picked case.
