#!/bin/bash -l

echoerr() { printf "%s\n" "$*" >&2; }

x1="0"
if [ -d "$1" ]; then
if [ -s "$1""/EXTRACT_1A-S2_target.fasta" ]; then
if [ -s "$1""/EXTRACT_2B-S2_target.fasta" ]; then
if [ -s "$1""/EXTRACT_2X-S2_target.fasta" ]; then
x1="1"
fi
fi
fi
fi

if [ "$x1" == "1" ]; then
  path=$2
  d1=$1
  echoerr "Data folder is $d1"
else
  echoerr "usage bash ./PBP_AA_sampledir_to_MIC_20171207 data_dir"
  echoerr ""
  echoerr "data_dir is a directory that must conatin 3 files with the following exact names, respectively:"
  echoerr "EXTRACT_1A-S2_target.fasta"
  echoerr "EXTRACT_2X-S2_target.fasta"
  echoerr "EXTRACT_2B-S2_target.fasta"
  echoerr ""
  echoerr "See README.txt for details"
  echoerr "Program not run"
  exit 1
fi

data_dir=${path}/bLactam_MIC_Rscripts

AAseqDir="$d1""/PBP_to_MIC_temp"
mkdir -p "${AAseqDir}"
cd  "${AAseqDir}" || exit

rm -f temp*
cp -f "$d1""/EXTRACT_1A-S2_target.fasta" temp1.fna
transeq.pl temp1.fna temp1.faa -frame=1
echo ">Sample1" > Sample_PBP1A_AA.faa
grep -v ">" temp1.faa >> Sample_PBP1A_AA.faa


rm -f temp*
cp -f "$d1""/EXTRACT_2B-S2_target.fasta" temp1.fna
transeq.pl temp1.fna temp1.faa -frame=1
echo ">Sample1" > Sample_PBP2B_AA.faa
grep -v ">" temp1.faa >> Sample_PBP2B_AA.faa

rm -f temp*
cp -f "$d1""/EXTRACT_2X-S2_target.fasta" temp1.fna
transeq.pl temp1.fna temp1.faa -frame=1
echo ">Sample1" > Sample_PBP2X_AA.faa
grep -v ">" temp1.faa >> Sample_PBP2X_AA.faa

rm -f temp*

#
scr1="${data_dir}/AAtoMICwrapper_2.sh"
bash "${scr1}" "${AAseqDir}" "${data_dir}"

#
fin="$AAseqDir"/Sample_PBPtype_MIC2_Prediction.csv
scr1="$path/bLactam_MIC_Rscripts/MIC_format_with_SIR.R"
Rscript "${scr1}" "${fin}"
fout="$AAseqDir"/Sample_PBPtype_MIC2_Prediction.csv_MIC_formatted_with_SIR.csv

#20180710
#Add NAs for 4 more beta-lactams
#WGS_AMP_SIGN WGS_AMP WGS_AMP_SIR 
#WGS_CPT_SIGN WGS_CPT WGS_CPT_SIR 
#WGS_ZOX_SIGN WGS_ZOX WGS_ZOX_SIR 
#WGS_FOX_SIGN WGS_FOX WGS_FOX_SIR 

cd "$AAseqDir" || exit
awk 'NR<=2'  $fout | awk -F"," '{$1=""; print $0}'  | sed 's/"//g' | sed 's/ //' > temp1.txt
echo "WGS_AMP_SIGN WGS_AMP WGS_AMP_SIR \
WGS_CPT_SIGN WGS_CPT WGS_CPT_SIR \
WGS_ZOX_SIGN WGS_ZOX WGS_ZOX_SIR \
WGS_FOX_SIGN WGS_FOX WGS_FOX_SIR" > temp2.txt

echo "NA NA NA \
NA NA NA \
NA NA NA \
NA NA NA " >> temp2.txt
paste -d ' ' temp1.txt temp2.txt>"$d1""/BLACTAM_MIC_RF_with_SIR.txt"

echoerr "BLACTAM MIC output file:" "${d1}/BLACTAM_MIC_RF_with_SIR.txt"

rm -f temp* 
#rm -rf $AAseqDir

#module unload EMBOSS/6.4.0
