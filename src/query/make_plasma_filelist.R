# Need to make a list of plasma metQTL files, since METSIM and Schlosser
# metabolites do not overlap
library(dplyr)
library(tidyr)

dat = read.delim('data/metqtl_pheno/plasma_metqtl_pheno.txt', header = T, stringsAsFactors = F)

# Generate file names:
metsim = dat %>% filter(is.na(summaryStatistics)) %>% pull(phenocode)
schl = dat %>% filter(is.na(phenocode)) %>% pull(summaryStatistics)
both = dat %>% filter(!(is.na(phenocode) | is.na(summaryStatistics))) %>% mutate(name = paste(phenocode, summaryStatistics, sep = '_')) %>% pull(name)

metsim = paste('data/metsim_data/', metsim, '.clean.txt.gz', sep = '')
schl = paste('data/schlosser_metqtl/plasma/', schl, '_buildGRCh37.tsv.gz', sep = '')
both = paste('data/plasma_meta/', both, '.meta1.tbl.gz', sep = '')

# Save list
writeLines(c(metsim, schl, both), 'data/metqtl_pheno/all_plasma_files.txt')
