#!/usr/bin/env bash
set -euo pipefail

script_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PATH="${script_root}:${PATH}"

assembly=/tmp/sequence.fa
allDB_dir="${script_root}/SPN_Reference_DB"
sample_out=$(mktemp -d /tmp/spn-pbp.XXXXXX)
trap 'rm -rf "$sample_out"' EXIT
temp_path="${script_root}"
cd "$sample_out"

echoerr() { printf "%s\n" "$*" >&2; }

###Start Doing Stuff###
#mkdir -p ${sample_out}
#cd "$sample_out"
just_name=$(basename "$sample_out")

###Call GBS bLactam Resistances###
PBP-Gene_Typer.pl -f "$assembly" -r "${allDB_dir}/MOD_bLactam_resistance.fasta" -o "$sample_out" -n "$just_name" -s SPN -p 1A,2B,2X

###Predict bLactam MIC###
scr1="${temp_path}/bLactam_MIC_Rscripts/PBP_AA_sampledir_to_MIC_20180710.sh"
bash "${scr1}" "$sample_out" "$temp_path"

if [[ -n "${SPN_PBP_DEBUG_DIR:-}" ]]; then
    debug_name="${SPN_PBP_DEBUG_NAME:-$just_name}"
    debug_out="${SPN_PBP_DEBUG_DIR}/${debug_name}"
    mkdir -p "$debug_out"
    cp EXTRACT_*.fasta "$debug_out/"
    cp PBP_to_MIC_temp/Sample_PBP*_AA.faa "$debug_out/"
    cp TEMP_pbpID_Results.txt BLACTAM_MIC_RF_with_SIR.txt "$debug_out/"
fi

###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
tabl_out="TABLE_Isolate_Typing_results.txt"
printf "${just_name}\t" >> "${tabl_out}"

###PBP_ID Output###
justPBPs="NF"
sed 1d TEMP_pbpID_Results.txt | while read -r line
do
    if [[ -n "$line" ]]
    then
        justPBPs=$(echo "$line" | awk -F"\t" '{print $2}' | tr ':' '\t')
    fi
    printf "%s\t" "$justPBPs" >> "$tabl_out"
done

pbpID=$(tail -n1 "TEMP_pbpID_Results.txt" | awk -F"\t" '{print $2}')
if [[ ! "$pbpID" =~ .*NF.* ]] #&& [[ ! "$pbpID" =~ .*NEW.* ]]
then
    echoerr "No NF outputs for PBP Type"
    bLacTab=$(tail -n1 "BLACTAM_MIC_RF_with_SIR.txt" | tr ' ' '\t')
    printf "%s\t" "$bLacTab" >> "$tabl_out"
else
    echoerr "One of the PBP types has an NF"
    printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\t" >> "$tabl_out"
fi

cat ${tabl_out}
