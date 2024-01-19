# Need to figure out which summary stats from METSIM and Schlosesr match in
# terms of the metabolite it is measuring
library(dplyr)
library(tidyr)

# Load METSIM
metsim = read.delim('data/metsim_data/phenotypes.tsv', header = T, stringsAsFactors = F)
metsim = metsim %>% select(phenocode:phenostring) %>% mutate(phenostring = gsub(' - ', '-', phenostring))

# Load in Schlosser
schl = read.table('data/schlosser_metqtl/PMID37277652_studies_export.tsv', sep = '\t', header = T, stringsAsFactors = F)

plasma = schl %>% filter(grepl('^Plasma', reportedTrait)) %>% select(reportedTrait:efoTraits, discoverySampleAncestry, associationCount, summaryStatistics)
urine = schl %>% filter(grepl('^Urine', reportedTrait)) %>% select(reportedTrait:efoTraits, discoverySampleAncestry, associationCount, summaryStatistics)

# For the Schlosser plasma phenotypes, some are just "metabolite measurement"
# in efiTraits, when there is a compound name in "reportedTrait" - pull the
# actual compound out:
plasma = plasma %>% mutate(phenostring = gsub(' levels in chronic kidney disease', '', reportedTrait), phenostring = gsub('^Plasma ', '', phenostring))
urine = urine %>% mutate(phenostring = gsub(' levels in chronic kidney disease', '', reportedTrait), phenostring = gsub('^Urine ', '', phenostring))

# What's the overlap:
which(metsim$phenostring %in% plasma$phenostring) %>% length # 1101 metabolites

# Make a full table:
full_table = plasma %>% select(phenostring, summaryStatistics) %>% mutate(summaryStatistics = gsub('.*/', '', summaryStatistics)) %>% full_join(metsim, .) %>% select(phenostring, phenocode, summaryStatistics)

# Save the full phenotype table:
write.table(full_table, 'data/metqtl_pheno/plasma_metqtl_pheno.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Just those that overlapped:
overlap = full_table %>% filter(!(is.na(phenocode) | is.na(summaryStatistics)))

# Save the overlap list so it can be used for meta-analysing
write.table(overlap, 'data/plasma_meta/plasma_metqtl_meta_list.txt', sep = '\t', col.names = F, row.names = F, quote = F)

# Save Schlosser phenotype table:
plasma = plasma %>% mutate(summaryStatistics = gsub('.*/', '', summaryStatistics))
urine = urine %>% mutate(summaryStatistics = gsub('.*/', '', summaryStatistics))

write.table(plasma, 'data/metqtl_pheno/schlosser_plasma_pheno.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(urine, 'data/metqtl_pheno/schlosser_urine_pheno.txt', sep = '\t', col.names = T, row.names = F, quote = F)
