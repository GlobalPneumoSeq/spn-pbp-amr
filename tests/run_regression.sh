#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
image=${IMAGE:-spn_pbp_amr:regression}
work=$(mktemp -d /tmp/spn-pbp-regression.XXXXXX)
cleanup() {
    rm -rf "$work"
}
trap cleanup EXIT

docker build --rm -t "$image" "$root" >/dev/null

Rscript "$root/tests/test_blosum_utils.R"
perl "$root/tests/test_extract_gene.pl"

if printf '' | docker run --rm -i --user "$(id -u):$(id -g)" "$image" >"$work/empty.out" 2>"$work/empty.err"; then
    echo "empty FASTA unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'ERROR: input must be a non-empty FASTA' "$work/empty.err"

if printf 'not fasta\n' | docker run --rm -i --user "$(id -u):$(id -g)" "$image" >"$work/malformed.out" 2>"$work/malformed.err"; then
    echo "malformed FASTA unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'ERROR: input must be a non-empty FASTA' "$work/malformed.err"

declare -A expected_status=(
    [SAMN20398278]='2:6:0'
    [SAMN20409264]='0:1:1'
    [SAMN33825264]='0:1:1'
)

run_fixture() {
    local fasta=$1 name actual fixture_work
    name=$(basename "$fasta" .fasta)
    actual="$work/${name}.json"
    fixture_work="$work/$name"
    mkdir -p "$fixture_work"

    docker run --rm -i --user "$(id -u):$(id -g)" \
        -e SPN_PBP_WORK_DIR="/work/$name" \
        -v "$work:/work" \
        "$image" < "$fasta" >"$actual" 2>"$work/${name}.err"
}

pids=()
for fasta in "$root"/tests/fixtures/examples/*.fasta; do
    run_fixture "$fasta" &
    pids+=("$!")
done
for pid in "${pids[@]}"; do
    wait "$pid"
done

for fasta in "$root"/tests/fixtures/examples/*.fasta; do
    name=$(basename "$fasta" .fasta)
    actual="$work/${name}.json"
    fixture_work="$work/$name"

    jq empty "$actual"
    diff -u \
        <(jq -S . "$actual") \
        <(jq -S . "$root/tests/fixtures/examples/expected-${name}.json")

    for artifact in EXTRACT_1A-S2_target.fasta EXTRACT_2B-S2_target.fasta EXTRACT_2X-S2_target.fasta \
        Sample_PBP1A_AA.faa Sample_PBP2B_AA.faa Sample_PBP2X_AA.faa; do
        if [[ "$artifact" == Sample_* ]]; then
            actual_artifact="$fixture_work/PBP_to_MIC_temp/$artifact"
        else
            actual_artifact="$fixture_work/$artifact"
        fi
        expected_artifact="$root/tests/fixtures/examples/expected-artifacts/$name/$artifact"
        if [[ "$artifact" == EXTRACT_* ]]; then
            # Historical reverse-strand extraction omitted the final newline,
            # whereas some baseline artifacts retained it. Compare the FASTA
            # content, not that presentation-only difference.
            diff -u \
                <(perl -0777 -pe 's/\n\z//' "$actual_artifact") \
                <(perl -0777 -pe 's/\n\z//' "$expected_artifact")
        else
            cmp "$actual_artifact" "$expected_artifact"
        fi
    done

    actual_status=$(awk 'NR == 2 { print $2 }' "$fixture_work/TEMP_pbpID_Results.txt")
    [[ "$actual_status" == "${expected_status[$name]}" ]]
done

run_fixture "$root/tests/fixtures/nf/no-pbp.fasta"
jq empty "$work/no-pbp.json"
diff -u \
    <(jq -S . "$work/no-pbp.json") \
    <(jq -S . "$root/tests/fixtures/nf/expected-no-pbp.json")
[[ "$(awk 'NR == 2 { print $2 }' "$work/no-pbp/TEMP_pbpID_Results.txt")" == "NF:NF:NF" ]]

# Derive a versioned novel-PBP fixture from a known assembly. The changed base
# is in PBP2B codon 1, so it must be reported as NEW while preserving a valid
# end-to-end JSON prediction.
novel_fasta="$work/novel-pbp.fasta"
perl -e '
    my ($contig, $position) = ("NODE_1_length_271402_cov_32.208704", 200991);
    my ($current, $offset, $changed);
    while (<>) {
        if (/^>(\S+)/) { $current = $1; $offset = 0; print; next; }
        if (defined $current && $current eq $contig) {
            my $bases = $_;
            $bases =~ s/\s+//g;
            my $length = length $bases;
            if ($offset <= $position && $position < $offset + $length) {
                my $column = $position - $offset;
                substr($_, $column, 1) = "A";
                $changed = 1;
            }
            $offset += $length;
        }
        print;
    }
    die "novel-PBP fixture mutation was not applied\n" unless $changed;
' "$root/tests/fixtures/examples/SAMN20398278.fasta" > "$novel_fasta"
run_fixture "$novel_fasta"
jq empty "$work/novel-pbp.json"
diff -u \
    <(jq -S . "$work/novel-pbp.json") \
    <(jq -S . "$root/tests/fixtures/novel/expected-novel-pbp.json")
[[ "$(awk 'NR == 2 { print $2 }' "$work/novel-pbp/TEMP_pbpID_Results.txt")" == "2:NEW:0" ]]

# Insert one base into PBP2B to shift its translation frame. Extraction still
# succeeds, but allele typing must report ERROR and the predictor must return
# a valid JSON result with that status.
error_fasta="$work/error-pbp.fasta"
perl -e '
    my ($contig, $position) = ("NODE_1_length_271402_cov_32.208704", 200995);
    my ($current, $offset, $changed);
    while (<>) {
        if (/^>(\S+)/) { $current = $1; $offset = 0; print; next; }
        if (defined $current && $current eq $contig) {
            my $bases = $_;
            $bases =~ s/\s+//g;
            my $length = length $bases;
            if (!$changed && $offset <= $position && $position < $offset + $length) {
                my $column = $position - $offset;
                substr($_, $column, 0) = "A";
                $changed = 1;
            }
            $offset += $length;
        }
        print;
    }
    die "error-PBP fixture mutation was not applied\n" unless $changed;
' "$root/tests/fixtures/examples/SAMN20398278.fasta" > "$error_fasta"
run_fixture "$error_fasta"
jq empty "$work/error-pbp.json"
diff -u \
    <(jq -S . "$work/error-pbp.json") \
    <(jq -S . "$root/tests/fixtures/error/expected-error-pbp.json")
[[ "$(awk 'NR == 2 { print $2 }' "$work/error-pbp/TEMP_pbpID_Results.txt")" == "2:ERROR:0" ]]

version=$(docker image inspect --format '{{index .Config.Labels "org.opencontainers.image.version"}}' "$image")
[[ "$version" == "0.2.0" ]]

docker run --rm --entrypoint /bin/sh "$image" -c '
    ! command -v bedtools
    ! command -v needle
    ! command -v curl
    Rscript -e "if (requireNamespace(\"glmnet\", quietly=TRUE)) quit(status=1)"
    Rscript -e "stopifnot(as.character(packageVersion(\"randomForest\")) == \"4.6.14\")"
'

echo "Docker regression tests passed"
