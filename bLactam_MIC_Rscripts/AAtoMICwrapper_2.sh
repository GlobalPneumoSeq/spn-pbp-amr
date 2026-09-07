#!/bin/bash

echoerr() { printf "%s\n" "$*" >&2; }

x1="0"
if [ -d "$1" ]; then
if [ -s "$1""/Sample_PBP1A_AA.faa" ]; then
if [ -s "$1""/Sample_PBP2B_AA.faa" ]; then
if [ -s "$1""/Sample_PBP2X_AA.faa" ]; then
x1="1"
fi
fi
fi
fi

if [ "$x1" == "1" ]; then
  AAseqDir=$1
  data_dir=$2
  echoerr "Data folder is $AAseqDir"
else
  echoerr "usage bash ./AAtoMICwrapper.sh data_dir"
  echoerr ""
  echoerr "data_dir is a directory that must conatin 3 files with the following exact names, respectively:"
  echoerr "Sample_PBP1A_AA.faa"
  echoerr "Sample_PBP2B_AA.faa"
  echoerr "Sample_PBP2X_AA.faa"
  echoerr ""
  echoerr "See README.txt for details"
  echoerr "Program not run"
  exit 1
fi

export SPN_PBP_DATA_DIR="$data_dir"
export SPN_PBP_MODEL_DIR="${data_dir}/newDB"

#
faaDir=${AAseqDir}"/Sample_AAtoMIC/faa/"
rm -rf "${faaDir}"
mkdir -p "${faaDir}"
cd "${faaDir}" || exit

cp "${data_dir}"/Ref_PBP_3.faa .
cp "${AAseqDir}"/*.faa .

#module load clustal-omega/1.2
scr1="${data_dir}/Build_PBP_AA_tableR3.2.2.R"
Rscript "${scr1}" "${faaDir}"

predir=${AAseqDir}"/Sample_AAtoMIC/pre/"
rm -rf   "${predir}"
mkdir -p "${predir}"
cp ./Sample_PBP_AA_table.csv "${predir}"

cd "${predir}" || exit

scr1="${data_dir}/AAtable_To_MIC_MM_RF_EN_2.R"
Rscript "${scr1}" "${predir}"

cp Sample_PBPtype_MIC2_Prediction.csv "${AAseqDir}"


echoerr "MIC pridiction results are in file:"
echoerr "$AAseqDir""/Sample_PBPtype_MIC2_Prediction.csv"

