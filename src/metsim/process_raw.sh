#! /bin/bash

# Process the downloaded METSIM summary stats

zcat ${1} | cut -f1-10 | awk 'NR == 1 {print $0}; NR > 1 && $10 >= 0.001 && length($3) == 1 && length($4) == 1 {print $0}' > ${1%%.gz}
