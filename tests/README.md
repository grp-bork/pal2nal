# Test corpus

Golden outputs captured from the original implementation, vendored here
as `reference/pal2nal.v14.pl`, running under Perl 5.34. They are the
reference the Python port must reproduce.

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
omitted when empty) and `golden/<name>.status` (exit status).

The script echoes its input file names into the output header, so the
paths in `cases.tsv` are relative and the runner must be started from
this directory. 127 cases.

## What is covered

| Group | Cases |
|---|---|
| Output formats | `clustal`, `paml`, `fasta`, `codon`, and the default |
| Filters | `-nogap`, `-nomismatch`, `-blockonly` and their combinations |
| Web modes | `-html`, `-nostderr` and combinations |
| Input handling | CLUSTAL and FASTA peptide input, LF / CRLF / CR line endings, DNA split across several files, lower-case DNA, ID-matched and ID-matched-but-reordered DNA |
| Codon tables | all 17 supported tables, each probed with all 64 codons, plus cross-checks in both directions against the universal code |
| Biology | in-frame stop codons, all frame-shift notations from `howto.html`, UTRs and polyA tails, protein/DNA mismatch (warned and removed), clean baseline |
| Diagnostics | usage with no arguments, `-h`, invalid output format, invalid codon table, the `-output codon` conflicts, missing input file |
| Input validation | a gap typed as an en dash, a UTF-8 BOM, a latin-1 file, an embedded NUL, protein pasted into the DNA field, an unknown residue on its own and under `-nomismatch` and `-blockonly` |

## Codon table probes

`data/allcodons.nuc` holds all 64 codons in NCBI order. For each table N,
`data/allcodons_t<N>_pep.fasta` is that table's translation of those 64
codons, taken from the official NCBI `ncbieaa` strings rather than from
the Perl, so the probe is an independent reference.

Running table N's probe under `-codontable N` should convert silently.
16 of the 17 tables do.

### Known bug: `-codontable 10`

Table 10 is the Euplotid nuclear code, whose single deviation from the
universal code is **TGA = Cys**. `pal2nal.v14.pl` does not implement it:

    "C" => "((U|T)G(U|T|C|Y))"                    # TGT, TGC -- TGA missing
    "*" => "(((U|T)A(A|G|R))|((T|U)GA))"          # TGA still a stop

so `-codontable 10` behaves exactly like `-codontable 1`. This is
confirmed by two cases: `codontable_10_own_probe` warns on both TGA
columns, and `codontable_t1_under_10` is silent where it should have
complained.

The goldens record the behaviour as it is. The port has to decide
whether to reproduce the bug or implement the table correctly; see
`PORTING.md`.
