# Need to calculate sample size for each metabolite, since it could be
# different depending on which data set (and/or meta'ed) it came from
library(dplyr)
library(tidyr)

pheno = read.delim('data/metqtl_pheno/plasma_metqtl_pheno.txt', header = T, stringsAsFactors = F)
plasma = read.delim('data/metqtl_pheno/schlosser_plasma_pheno.txt', header = T, stringsAsFactors = F)
urine = read.delim('data/metqtl_pheno/schlosser_urine_pheno.txt', header = T, stringsAsFactors = F)

# Merge urine and plasma into a single table
schl = plasma %>% separate(discoverySampleAncestry, into = c('N_schlosser', 'ancestry')) %>% select(phenostring, summaryStatistics, N_schlosser)
schl = urine %>% separate(discoverySampleAncestry, into = c('N_schlosser', 'ancestry')) %>% select(phenostring, summaryStatistics, N_schlosser) %>% full_join(schl, ., by = 'phenostring', suffix = c('.plasma', '.urine'))

# To the schlosser table, merge the METSIM info
metsim = pheno %>% filter(!is.na(phenocode)) %>% mutate(N = 6136) %>% select(phenostring, phenocode, N)

full_pheno = full_join(schl, metsim, by = 'phenostring')
full_pheno[is.na(full_pheno)] = 0

full_pheno$N_plasma = as.numeric(full_pheno$N_schlosser.plasma) + full_pheno$N
full_pheno$N_urine = as.numeric(full_pheno$N_schlosser.urine)

full_pheno = full_pheno %>% select(phenostring, phenocode, contains('summary'), N_plasma, N_urine)

write.table(full_pheno, 'data/metqtl_pheno/full_pheno.txt', sep = '\t', col.names = T, row.names = F, quote = F)

