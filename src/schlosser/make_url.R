library(dplyr)
library(tidyr)

dat = read.table('data/schlosser_metqtl/PMID37277652_studies_export.tsv', sep = '\t', header = T, stringsAsFactors = F)

plasma = dat %>% filter(grepl('^Plasma', reportedTrait))
urine = dat %>% filter(grepl('^Urine', reportedTrait))

plasma_url = paste(plasma$summaryStatistics, '/', plasma$accessionId, '_buildGRCh37.tsv.gz', sep = '')
urine_url = paste(urine$summaryStatistics, '/', urine$accessionId, '_buildGRCh37.tsv.gz', sep = '')

writeLines(plasma_url, 'data/schlosser_metqtl/plasma_urls.txt')
writeLines(urine_url, 'data/schlosser_metqtl/urine_urls.txt')
