#!/bin/sh
#
# Regenerate the golden outputs by running the original Perl implementation.
#
#   sh generate_golden.sh [path-to-pal2nal.pl]
#
# Defaults to the last deployed version, reference/pal2nal.v14.pl.
# Run from this directory: the recorded output contains the input paths, so
# the relative paths in cases.tsv are what keep the goldens reproducible.

set -e
PAL2NAL=${1:-reference/pal2nal.v14.pl}
cd "$(dirname "$0")"
rm -rf golden
mkdir -p golden

while IFS="$(printf '\t')" read -r name args; do
    [ -n "$name" ] || continue
    # shellcheck disable=SC2086
    if perl "$PAL2NAL" $args > "golden/$name.out" 2> "golden/$name.err"; then
        status=0
    else
        status=$?
    fi
    printf '%s\n' "$status" > "golden/$name.status"
    [ -s "golden/$name.err" ] || rm -f "golden/$name.err"
done < cases.tsv

echo "$(ls golden/*.out | wc -l | tr -d ' ') cases written to golden/"
