#!/usr/bin/env bash
set -euo pipefail

cat - > /tmp/sequence.fa

if ! perl -ne '
    $records++ if /^>/;
    if ($seen_header) { $bases += length($_) unless /^>/; }
    $seen_header = 1 if /^>/;
    END { exit(!($records && $bases)); }
' /tmp/sequence.fa; then
    echo "ERROR: input must be a non-empty FASTA file containing at least one sequence" >&2
    exit 2
fi

trap 'rm -f /tmp/sequence.fa' EXIT

if ! pw_wrapper.sh | to_json.pl; then
    echo "ERROR: PBP prediction failed" >&2
    exit 1
fi
