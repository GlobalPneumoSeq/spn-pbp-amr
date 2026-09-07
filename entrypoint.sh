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

wrapper_output=$(mktemp /tmp/spn-pbp-output.XXXXXX)
trap 'rm -f "$wrapper_output" /tmp/sequence.fa' EXIT

if ! pw_wrapper.sh >"$wrapper_output"; then
    echo "ERROR: PBP prediction failed" >&2
    exit 1
fi

result_line=$(awk 'NF >= 25 { line = $0 } END { if (line) print line }' "$wrapper_output")
if [[ -z "$result_line" ]]; then
    echo "ERROR: predictor did not produce a complete result row" >&2
    exit 1
fi

printf '%s\n' "$result_line" | to_json.pl
