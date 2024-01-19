#! /bin/bash

# Script to make LD file for all the variants used in metQTL
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load plink/plink1.9b6.10
module load bcftools
module load htslib

CHR=$(echo $1 | sed 's/23/X/g')
START=$(echo $2 | awk '{print $1 - 500000}')
END=$(echo $2 | awk '{print $1 + 500000}')
SNP=$3
KGP=/Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/EUR/EUR_chr${CHR}.no_relatives.no_indel.biallelic.vcf.gz

bcftools view --regions ${CHR}:${START}-${END} --output-type z --output-file results/ld/${SNP}.vcf.gz ${KGP}

plink --vcf results/ld/${SNP}.vcf.gz --allow-no-sex --snps-only --r2 --inter-chr --ld-snp ${SNP} --ld-window-r2 0 --out results/ld/${SNP}

