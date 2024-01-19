library(dplyr)

dat = read.delim('data/metsim_data/phenotypes.tsv', header = T, stringsAsFactors = F)

url = paste('https://pheweb.org/metsim-metab/download/', dat$phenocode, sep = '') %>% unique

writeLines(url, 'data/metsim_data/dl_urls.txt')
