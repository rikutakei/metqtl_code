#! /bin/bash

module load gatk

# Create initial list of variants
zcat data/metsim_data/C98.tsv.gz | cut -f1-4 > data/metsim_data/tmp_var_list.txt

# Pull out variants not listed in the initial list
parallel -j50 "zgrep -Fwvf data/metsim_data/tmp_var_list.txt {} | cut -f1-4 > {.}.list" ::: $(ls data/metsim_data/C*.tsv.gz)

cat <(tail -n+2 data/metsim_data/tmp_var_list.txt) data/metsim_data/*.list > data/metsim_data/metsim_var.txt

awk '{print "chr"$1, $2, $1"_"$2, $3, $4, ".", ".", "."}' data/metsim_data/metsim_var.txt > data/metsim_data/tmp && mv data/metsim_data/tmp data/metsim_data/metsim_var.txt

cat data/liftover/vcf_header.txt data/metsim_data/metsim_var.txt | tr ' ' '\t' > data/metsim_data/metsim_var.vcf

gatk LiftoverVcf -I=data/metsim_data/metsim_var.vcf -O=data/metsim_data/metsim_var.hg19.vcf -C=data/liftover/hg38ToHg19.over.chain.gz --REJECT=data/metsim_data/rejected_var.vcf -R=/Volumes/archive/merrimanlab/reference_files/FASTA/hg19/hg19.fa

rm data/metsim_data/*.list
