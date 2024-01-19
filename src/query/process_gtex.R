# Pull out lead variant +/-500kb from the GTEx summary stats
library(vroom)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly = T)

gtex = read.table(args[1], header = T, sep = '\t', stringsAsFactors = F)
snp_list = vroom('data/eqtl/gtex_lookup_sorted_rsid.txt', col_names = c('gtex_id', 'rsid'))

gtex = left_join(gtex, snp_list, by = c('variant_id' = 'gtex_id'))

write.table(gtex, gsub('.tmp', '', args[1]), sep = '\t', col.names = T, row.names = F, quote = F)
