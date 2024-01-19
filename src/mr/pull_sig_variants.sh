#! /bin/bash
# Pull out most significant variant (i.e. lead SNP) from a summary stats

module load plink/plink1.9b6.10

plink1.9b6.10 --bfile data/1kgp_ref/${ANCESTRY}_wgs --clump ${OUT}.clump.in --clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0 --clump-kb 500 --out ${OUT}


