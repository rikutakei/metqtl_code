#! /bin/bash

# Pull out all of the relevant variants that are common across all of the data
# sets and format it to be ready for the MR analysis.
# Also prepare/generate summary for PLINK clumping

module load parallel
module load R

# Start with Major gout GWAS:
grep -Fwf data/mr_dat/sig_snps.common.txt data/gwas/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt > data/mr_dat/gout.sumstats.tmp
awk '{print $3, $1, $2, $4, $5, $6, $10, $11, $12}' data/mr_dat/gout.sumstats.tmp > data/mr_dat/gout.sumstats.txt

# METSIM:
parallel "grep -Fwf data/mr_dat/sig_snps.common.txt data/metsim_data/{}.tsv > data/mr_dat/metsim_{}.sumstats.tmp" ::: $(cat results/mr_results/metsim.mr_list.txt)
parallel "awk '{print \$5, \$1, \$2, \$3, \$4, \$10, \$8, \$9, \$7}' data/mr_dat/metsim_{}.sumstats.tmp > data/mr_dat/metsim_{}.sumstats.txt" ::: $(cat results/mr_results/metsim.mr_list.txt)

# Schlosser:
parallel "zgrep -Fwf data/mr_dat/sig_snps.common.txt data/schlosser_metqtl/plasma/{}_buildGRCh37.tsv.gz > data/mr_dat/schlosser_{}.sumstats.tmp" ::: $(cat results/mr_results/schl.mr_list.txt)
parallel "awk '{print \$1, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$2}' data/mr_dat/schlosser_{}.sumstats.tmp > data/mr_dat/schlosser_{}.sumstats.txt" ::: $(cat results/mr_results/schl.mr_list.txt)

# Chen:
parallel "grep -Fwf data/mr_dat/sig_snps.common.txt data/chen_metqtl/{}_buildGRCh38.tsv > data/mr_dat/chen_{}.sumstats.txt" ::: $(cat results/mr_results/chen.mr_list.txt)

rm data/mr_dat/*.tmp

# Add new header, since they are all in the same order
echo SNP CHR POS effect_allele other_allele MAF effect SE P > data/mr_dat/header.txt

parallel "cat data/mr_dat/header.txt {} | tr ' ' '\t' > {}.tmp && mv {}.tmp {}" ::: $(ls data/mr_dat/*sumstats.txt)
parallel "Rscript src/mr/prep_clump.R {}" ::: $(ls data/mr_dat/*sumstats.txt)

