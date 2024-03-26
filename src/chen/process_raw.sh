#! /bin/bash

# Process the downloaded METSIM summary stats

echo "SNP CHR POS effect_allele other_allele MAF effect SE P" | tr -s ' ' '\t' > ${1%%.gz}

zcat ${1} | awk 'NR > 1 && length($3) == 1 && length($4) == 1 {print $9, $0}' | tr -s ' ' '\t' | cut -f1-9 >> ${1%%.gz}
