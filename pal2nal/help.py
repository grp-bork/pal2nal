"""The usage text, captured from pal2nal.v14.pl:1942 rather than retyped.

The layout is v14's own, trailing space after "Show help" included, and the
only entries added since are --version (v15) and -partial (v16), which v14
has no equivalent of; leaving them out would mean shipping options that -h
denies exist. Everything else is
byte for byte the Perl's, so `tests/golden/help_flag.err` stays a usable
diff against it.

In v14 this always went to stderr even under -html, and -h did not stop the
run; v15 still writes to stderr but exits (see PORTING.md).
"""

from __future__ import annotations

import sys
from typing import TextIO

USAGE = """
pal2nal.pl  (v14)

Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]


    pep.aln:    protein alignment either in CLUSTAL or FASTA format

    nuc.fasta:  DNA sequences (single multi-fasta or separated files)

    Options:  -h            Show help 

              --version     Show the version and exit

              -output (clustal|paml|fasta|codon)
                            Output format; default = clustal

              -blockonly    Show only user specified blocks
                            '#' under CLUSTAL alignment (see example)

              -nogap        remove columns with gaps and inframe stop codons

              -nomismatch   remove mismatched codons (mismatch between
                            pep and cDNA) from the output

              -partial      when a peptide and its DNA do not correspond,
                            convert the codons that can be matched, gap the
                            rest, and carry on instead of stopping the run

              -codontable  N
                    1  Universal code (default)
                    2  Vertebrate mitochondrial code
                    3  Yeast mitochondrial code
                    4  Mold, Protozoan, and Coelenterate Mitochondrial code
                       and Mycoplasma/Spiroplasma code
                    5  Invertebrate mitochondrial code
                    6  Ciliate, Dasycladacean and Hexamita nuclear code
                    9  Echinoderm and Flatworm mitochondrial code
                   10  Euplotid nuclear code
                   11  Bacterial, archaeal and plant plastid code
                   12  Alternative yeast nuclear code
                   13  Ascidian mitochondrial code
                   14  Alternative flatworm mitochondrial code
                   15  Blepharisma nuclear code
                   16  Chlorophycean mitochondrial code
                   21  Trematode mitochondrial code
                   22  Scenedesmus obliquus mitochondrial code
                   23  Thraustochytrium mitochondrial code


              -html         HTML output (only for the web server)

              -nostderr     No STDERR messages (only for the web server)


    - sequence order in pep.aln and nuc.fasta should be the same.

    - IDs in pep.aln are used in the output.

"""


def showhelp(stream: TextIO | None = None) -> None:
    (stream if stream is not None else sys.stderr).write(USAGE)
