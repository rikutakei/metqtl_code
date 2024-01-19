#! /bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load plink/plink1.9b6.10
module load bcftools
module load htslib

plink --bfile data/1kgp_eur/EUR_wgs --maf 0.001 --allow-no-sex --r2 --ld-snp-list data/ld_proxy/all_sig_variants.txt --ld-window-r2 0.8 --ld-window 1000 --ld-window-kb 1000 --out results/proxy/1kgp_proxy

