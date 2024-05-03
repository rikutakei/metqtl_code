#! /bin/bash

module load R

# Make a list of independent lead variants for each metabolite
Rscript src/mr/select_mr_variants.R

module load parallel

# Make a list of all the variants and pull these out from the summary stats:
cat /scratch/rtakei/projects/2023/metqtl_project/data/mr_dat/*mr_snps.txt | sort | uniq > /scratch/rtakei/projects/2023/metqtl_project/data/mr_dat/all_mr_snps.txt

parallel "grep -Fwf /scratch/rtakei/projects/2023/metqtl_project/data/mr_dat/all_mr_snps.txt {} > {.}.mr.txt " ::: $(ls /scratch/rtakei/projects/2023/metqtl_project/data/mr_dat/*sumstats.txt)

