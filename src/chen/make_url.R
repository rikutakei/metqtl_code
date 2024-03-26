library(dplyr)
library(tidyr)

dat = read.table('data/chen_metqtl/PMID36635386_studies_export.tsv', sep = '\t', header = T, stringsAsFactors = F)

# Pull out EUR data:
dat = dat %>% filter(grepl('Eur', discoverySampleAncestry))

# Pull out the sample size for each trait
dat = dat %>% separate(discoverySampleAncestry, into = c('N', 'ancestry'), sep = ' ', remove = T, convert = T)

url = paste(dat$summaryStatistics, '/', dat$accessionId, '_buildGRCh38.tsv.gz', sep = '')
url = gsub('http', 'https', url)

writeLines(url, 'data/chen_metqtl/chen_urls.txt')
