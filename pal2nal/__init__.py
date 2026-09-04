"""PAL2NAL: convert a protein sequence alignment and the corresponding DNA
sequences into a codon alignment.

A Python port of pal2nal.pl v14 by Mikita Suyama. See CHANGELOG.md for the
changes made since v14 and PORTING.md for the decisions behind them.
"""

#: The single source of the version. pyproject.toml reads this at build
#: time (``[tool.setuptools.dynamic]``), so this is also what
#: ``importlib.metadata.version("pal2nal")`` reports for an installed copy.
__version__ = "16.0"
