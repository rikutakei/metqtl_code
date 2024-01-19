#! /bin/bash

# Script to concatenate the split 1000GP VCFs and convert it into PLINK bfiles

module load plink/plink1.9b6.10
module load bcftools/bcftools-1.11
module load htslib/htslib-1.11

DAT_PATH=/Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/EUR/

# Convert VCF into PLINK bfiles and add variant IDs:
parallel -j 22 "plink --vcf ${DAT_PATH}/EUR_chr{}.no_relatives.no_indel.biallelic.vcf.gz --double-id --set-missing-var-ids @_# --make-bed --out data/1kgp_eur/EUR_chr{}_tmp" ::: {1..22}

# One variant in X causes hassle, so remove it before you begin:
echo X 32362662 | tr ' ' '\t' > data/1kgp_eur/exclude_kgp_x.list

bcftools view -T ^data/1kgp_eur/exclude_kgp_x.list ${DAT_PATH}/EUR_chrX.no_relatives.no_indel.biallelic.vcf.gz -Oz -o data/1kgp_eur/EUR_chrX.no_relatives.no_indel.biallelic.clean.vcf.gz --threads 2

plink --vcf data/1kgp_eur/EUR_chrX.*.clean.vcf.gz --double-id --set-missing-var-ids @_#_\\\$1_\\\$2 --make-bed --out data/1kgp_eur/EUR_chrX_tmp

awk '/rs/ {print}; !/rs/ {split($2, chrpos, "_"); $2 = chrpos[1]"_"chrpos[2]; print}' data/1kgp_eur/EUR_chrX_tmp.bim | tr " " "\t" > data/1kgp_eur/tmp_x.bim && mv data/1kgp_eur/tmp_x.bim data/1kgp_eur/EUR_chrX_tmp.bim

# Remove any SNPs that are duplicated (all have same set of variants, so just
# use bim file from one ancestry):
cut -f2 data/1kgp_eur/EUR_chr{{1..22},X}_tmp.bim | sort | uniq -d > data/1kgp_eur/all_dup.txt

parallel -j 22 "plink --bfile data/1kgp_eur/EUR_chr{}_tmp --exclude data/1kgp_eur/all_dup.txt --make-bed --out data/1kgp_eur/EUR_chr{}" ::: {{1..22},X}

# Generate list of bfiles to merge:
ls data/1kgp_eur/EUR_chr{{1..22},X}.bim | sed 's/.bim//g' | sort -V > data/1kgp_eur/EUR_file_list.txt

plink --merge-list data/1kgp_eur/EUR_file_list.txt --out data/1kgp_eur/EUR_wgs

rm data/1kgp_eur/*tmp*
rm data/1kgp_eur/*nosex
