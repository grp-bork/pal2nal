"""Finding a peptide's codons in a DNA sequence.

pal2nal.v14.pl built one regex per aligned residue, concatenated them into
a single pattern and handed the result to Perl's engine (sub pn2codon,
lines 1159-1919). ``convert.py`` still works from exactly those patterns --
they are the transcribed genetic codes and stay the readable specification
-- but it no longer asks a regex engine to run them, for two reasons that
only show up at scale:

* the pattern is proportional to the peptide, so a 3000-residue alignment
  builds a 40 kB pattern per sequence and Python's ``re`` -- whose parser
  is written in Python -- spends far longer compiling it than matching it;
* an unanchored search is quadratic when the peptide nearly matches. A
  thousand alanines against a poly-GCT sequence re-tries a thousand codons
  at every offset before failing.

Both go away once you notice what the patterns actually are. **Every
fragment matches exactly three characters** -- there is no quantifier, no
anchor and no backreference anywhere in the thirty-one tables, and
``tests/test_matcher.py`` fails if a table added later breaks that. So the
concatenated pattern is fixed-width, "the leftmost match" is simply the
smallest offset at which every fragment matches its own codon, and
backtracking cannot change the answer.

That turns the search into a bit-parallel AND. One integer per alphabet
character records the positions that character occupies in the DNA, built
with ``str.translate`` so the per-nucleotide work happens in C. A
fragment's positions are ORs and ANDs of those, computed once per distinct
fragment. Matching a whole peptide is then one shifted AND per residue,
on integers, stopping the moment no offset survives.

The bit indices are DNA positions throughout: bit *i* of an offset set
means "this matches starting at nucleotide *i*". Shifting a set right by
*d* therefore asks "does it match *d* further along", which is why the
whole search is a loop of ``acc &= positions >> offset``.
"""

from __future__ import annotations

from functools import cache
from typing import Final

#: One position of one codon: the characters allowed there, or None for
#: ``.``, which matches anything at all. None is not the same as "every
#: base": a degenerate third position matches an IUPAC code no table
#: mentions, and v14's regexes did too, so the distinction is load-bearing.
#: A parsed fragment is a tuple of the alternative codons it accepts, each
#: of them exactly three of these.
Atom = frozenset[str] | None
Codon = tuple[Atom, Atom, Atom]


def _split(text: str, sep: str = "|") -> list[str]:
    """Split on `sep` at paren depth zero."""
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    parts.append("".join(current))
    return parts


def _group(text: str, start: int) -> int:
    """Index just past the ")" closing the "(" at `start`."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return i + 1
    raise ValueError(f"unbalanced group in codon pattern: {text!r}")


def _parse(text: str) -> list[list[Atom]]:
    """One codon pattern into alternatives of Atoms.

    The grammar is only groups, ``|``, literal bases and ``.``, so this is
    a dozen lines rather than a regex engine. Concatenation multiplies out:
    a group nested inside an alternative contributes each of its own
    alternatives. A group whose alternatives are all a single character is
    folded into one Atom instead -- that is the ``(A|G|R)`` of the tables
    read as the character class it really is, and it keeps the number of
    alternatives at one or two for almost every fragment.
    """
    out: list[list[Atom]] = []
    for alternative in _split(text):
        seqs: list[list[Atom]] = [[]]
        i = 0
        while i < len(alternative):
            ch = alternative[i]
            if ch == "(":
                end = _group(alternative, i)
                inner = _parse(alternative[i + 1 : end - 1])
                if all(len(s) == 1 for s in inner):
                    # a character class written as an alternation
                    merged: set[str] = set()
                    degenerate = False
                    for s in inner:
                        if s[0] is None:
                            degenerate = True
                        else:
                            merged |= s[0]
                    atom = None if degenerate else frozenset(merged)
                    seqs = [s + [atom] for s in seqs]
                else:
                    seqs = [s + extra for s in seqs for extra in inner]
                i = end
            else:
                atom = None if ch == "." else frozenset(ch)
                seqs = [s + [atom] for s in seqs]
                i += 1
        out.extend(seqs)
    return out


@cache
def branches(fragment: str) -> tuple[Codon, ...]:
    """`_parse`, memoised: a table has at most two dozen fragments and a
    run touches them over and over."""
    parsed = _parse(fragment)
    if any(len(branch) != 3 for branch in parsed):
        raise ValueError(f"codon pattern is not three wide: {fragment!r}")
    return tuple(tuple(branch) for branch in parsed)


def codon_matches(fragment: str, codon: str) -> bool:
    """Whether one three-character codon is one this fragment accepts.

    Case-insensitive, like the ``re.IGNORECASE`` search it replaces, and
    false for anything that is not three characters -- which is what the
    old ``re.search`` returned for the short slice at a sequence's end.
    """
    if len(codon) != 3:
        return False
    upper = codon.upper()
    return any(
        all(atom is None or upper[i] in atom for i, atom in enumerate(branch))
        for branch in branches(fragment)
    )


class Offsets:
    """Where each codon fragment matches in one DNA sequence.

    Construction is one pass per distinct character; everything after that
    is integer arithmetic, so the cost of a search is the peptide's length
    times the sequence's length in *words*, not in characters.
    """

    __slots__ = ("_any", "_at", "_cache", "length")

    def __init__(self, nuc: str) -> None:
        self.length = len(nuc)
        upper = nuc.upper()
        present = set(upper)
        #: positions of each character, as a bit per nucleotide. The
        #: reversal puts nucleotide 0 in bit 0, so ">> d" means "d further
        #: along the sequence".
        self._at: dict[str, int] = {}
        for ch in present:
            marks = {ord(c): ("1" if c == ch else "0") for c in present}
            self._at[ch] = int(upper.translate(marks)[::-1], 2)
        #: every position; also what "." matches, since it matches any
        #: character and not merely any base
        self._any = (1 << self.length) - 1
        self._cache: dict[str, int] = {}

    def _atom(self, atom: Atom) -> int:
        if atom is None:
            return self._any
        positions = 0
        for ch in atom:
            positions |= self._at.get(ch, 0)
        return positions

    def fragment(self, pattern: str) -> int:
        """The offsets at which one codon fragment matches.

        A codon that would run off the end drops out for free: the third
        position's set is shifted down by two, so the top two bits clear.
        """
        cached = self._cache.get(pattern)
        if cached is not None:
            return cached
        found = 0
        for branch in branches(pattern):
            here = self._any
            for i, atom in enumerate(branch):
                here &= self._atom(atom) >> i
                if not here:
                    break
            found |= here
        self._cache[pattern] = found
        return found

    def wildcard(self, width: int) -> int:
        """Offsets where `width` characters of anything still fit: what a
        frame-shift numeral consumes."""
        return self._any >> (width - 1) if width else self._any

    def search(self, steps: list[tuple[int, str | None]]) -> int | None:
        """The leftmost offset at which every step matches in turn, or None.

        `steps` pairs a width with a fragment, or with None for a
        frame-shift numeral, which constrains nothing but still advances
        the offset. An empty `steps` is the caller's business, not this
        method's: it reports offset 0, matching nothing.
        """
        surviving = self._any
        offset = 0
        cache = self._cache          # one residue per step: keep it local
        for width, pattern in steps:
            if pattern is None:
                positions = self.wildcard(width)
            else:
                positions = cache.get(pattern)
                if positions is None:
                    positions = self.fragment(pattern)
            surviving &= positions >> offset
            if not surviving:
                return None
            offset += width
        # the lowest set bit is the leftmost surviving offset
        return (surviving & -surviving).bit_length() - 1


__all__: Final = ["Offsets", "branches", "codon_matches"]
