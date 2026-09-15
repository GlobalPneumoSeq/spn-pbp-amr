#!/usr/bin/env bash
set -euo pipefail

assembly=$1
allDB_dir="/predictor/SPN_Reference_DB"
if [[ -n "${SPN_PBP_WORK_DIR:-}" ]]; then
    sample_out=$SPN_PBP_WORK_DIR
    mkdir -p "$sample_out"
    cleanup_workdir=false
else
    sample_out=$(mktemp -d /tmp/spn-pbp.XXXXXX)
    cleanup_workdir=true
fi
cleanup() {
    if "$cleanup_workdir"; then
        rm -rf "$sample_out"
    fi
}
trap cleanup EXIT
cd "$sample_out"

echoerr() { printf "%s\n" "$*" >&2; }

###Start Doing Stuff###
#mkdir -p ${sample_out}
#cd "$sample_out"
just_name=$(basename "$sample_out")

###Call GBS bLactam Resistances###
PBP-Gene_Typer.pl -f "$assembly" -r "${allDB_dir}/MOD_bLactam_resistance.fasta" -o "$sample_out" -n "$just_name" -s SPN -p 1A,2B,2X >&2

###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
tabl_out="TABLE_Isolate_Typing_results.txt"

###PBP_ID Output###
pbpID=$(tail -n1 "TEMP_pbpID_Results.txt" | awk -F"\t" '{print $2}')
IFS=: read -r pbp1a pbp2b pbp2x <<< "$pbpID"
printf '%s\t%s\t%s\t%s\t' "predictor" "$pbp1a" "$pbp2b" "$pbp2x" >> "$tabl_out"

if [[ "$pbpID" != *NF* ]] #&& [[ ! "$pbpID" =~ .*NEW.* ]]
then
    echoerr "No NF outputs for PBP Type"
    ###Predict bLactam MIC###
    bash PBP_AA_sampledir_to_MIC_20180710.sh "$sample_out" "$sample_out" >&2
    bLacTab=$(tail -n1 "BLACTAM_MIC_RF_with_SIR.txt" | tr ' ' '\t')
    printf "%s\t" "$bLacTab" >> "$tabl_out"
else
    echoerr "One of the PBP types has an NF"
    printf 'NF\t%.0s' {1..21} >> "$tabl_out"
fi

cat "$tabl_out"
