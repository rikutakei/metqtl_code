#! /bin/bash

# Pull out relevant GTEx gene data and process it so it includes variants in
# the region of interest

GTEX_ID=$1
CHR=$(echo ${GTEX_ID} | sed 's/_.*//g')
ENSG=$2
TISSUE=$3
SNP=$4

cat <(head -1 /Volumes/archive/merrimanlab/reference_files/GTEx_v8_resources/GTEx_eQTL_v8/GTEx_files_by_chromosome/${TISSUE}_${CHR}.allpairs.txt) <(grep -Fw ${ENSG} /Volumes/archive/merrimanlab/reference_files/GTEx_v8_resources/GTEx_eQTL_v8/GTEx_files_by_chromosome/${TISSUE}_${CHR}.allpairs.txt) > results/eqtl/${SNP}_${TISSUE}_${ENSG%%.*}.tmp.txt

Rscript src/query/process_gtex.R results/eqtl/${SNP}_${TISSUE}_${ENSG%%.*}.tmp.txt

rm results/eqtl/${SNP}_${TISSUE}_${ENSG%%.*}.tmp.txt

