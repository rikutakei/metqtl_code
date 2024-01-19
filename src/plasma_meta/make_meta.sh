#! /bin/bash

METSIM=data/metsim_data/${1}.tsv.gz
SCHL=data/schlosser_metqtl/plasma/${2}_buildGRCh37.tsv.gz
OUT_NAME=data/plasma_meta/${1}_${2}.meta

# METSIM
Rscript src/plasma_meta/clean_metsim.R ${METSIM}

# Schlosser
cat <(echo cpid marker_id P CHR POS effect_allele other_allele MAF effect SE) <(zcat ${SCHL} | tail -n+2 | awk '{print $3"_"$4, $0}') | tr ' ' '\t' > ${SCHL%%.*}.metal.txt

# Generate the METAL file
sed -e "s,METSIM_SUMSTATS,${METSIM%%.*}.clean.txt,g" -e "s,SCHLOSSER_SUMSTATS,${SCHL%%.*}.metal.txt,g" -e "s,META_RESULT,${OUT_NAME},g" src/plasma_meta/metal_template.txt > data/plasma_meta/${1}_${2}.metal.txt

# Run METAL
metal data/plasma_meta/${1}_${2}.metal.txt > data/plasma_meta/${1}_${2}.metal.log

# METSIM data takes a while to generate (and also it will be in hg19, so it
# would be useful to keep for later), so just remove the temporary Schlosesr
# file and gzip the METSIM and meta-analysis files
rm ${SCHL%%.*}.metal.txt
rm ${METSIM%%.*}.clean.txt.gz
gzip ${METSIM%%.*}.clean.txt
gzip ${OUT_NAME}1.tbl
