#! /bin/bash

module load parallel

# Pull out all the relevant rsIDs from each of the summary stats and find all
# of the common variants.
# Then use these variants for the MR analysis

# First make a list of variants that are P <= 1e-3 in any of the summary stats

# Start with Major gout GWAS:
awk '$22 >= 3 {print $3}' data/gwas/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt > data/mr_dat/gout.sig_snps.txt

# METSIM:
parallel "awk '\$7 <= 0.001 {print \$5}' data/metsim_data/{}.tsv > data/mr_dat/metsim_{}.sig_snps.txt" ::: $(cat results/mr_results/metsim.mr_list.txt)

# Schlosser:
parallel "zcat data/schlosser_metqtl/plasma/{}_buildGRCh37.tsv.gz | awk '\$2 <= 0.001 {print \$1}' > data/mr_dat/schl_{}.sig_snps.txt" ::: $(cat results/mr_results/schl.mr_list.txt)

# Chen:
parallel "awk '\$9 <= 0.001 {print \$1}' data/chen_metqtl/{}_buildGRCh38.tsv > data/mr_dat/chen_{}.sig_snps.txt" ::: $(cat results/mr_results/chen.mr_list.txt)

cat data/mr_dat/*sig_snps.txt | grep rs | sort --parallel=2 | uniq > data/mr_dat/tmp && mv data/mr_dat/tmp data/mr_dat/sig_snps.txt


# Now figure out what variants are present in all of the summary stats:
grep -Fwf data/mr_dat/sig_snps.txt data/gwas/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt | cut -f3 | sort > data/mr_dat/gout.found.txt

# METSIM:
parallel "grep -Fwf data/mr_dat/sig_snps.txt data/metsim_data/{}.tsv | cut -f5 | sort > data/mr_dat/metsim_{}.found.txt" ::: $(cat results/mr_results/metsim.mr_list.txt)

# Schlosser:
parallel "zgrep -Fwf data/mr_dat/sig_snps.txt data/schlosser_metqtl/plasma/{}_buildGRCh37.tsv.gz | cut -f1 | sort > data/mr_dat/schl_{}.found.txt" ::: $(cat results/mr_results/schl.mr_list.txt)

# Chen:
parallel "grep -Fwf data/mr_dat/sig_snps.txt data/chen_metqtl/{}_buildGRCh38.tsv | cut -f1 | sort > data/mr_dat/chen_{}.found.txt" ::: $(cat results/mr_results/chen.mr_list.txt)

cp data/mr_dat/sig_snps.txt data/mr_dat/sig_snps.common.txt

for i in $(ls data/mr_dat/*.found.txt); do
comm -12 data/mr_dat/sig_snps.common.txt ${i} > data/mr_dat/tmp && mv data/mr_dat/tmp data/mr_dat/sig_snps.common.txt
done

