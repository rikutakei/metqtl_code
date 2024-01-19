library(dplyr)
library(tidyr)

dat = read.table('data/gwas/indep_snps_from_paper.txt', sep = '\t', header = T, stringsAsFactors = F, skip = 1)
dat = dat %>% filter(European == 'European') %>% select(X, rsID:GWAS.Sex)
colnames(dat) = c('locus', 'SNP', 'chrpos', 'cohort')
dat = dat %>% separate(chrpos, sep = ':', into = c('chr', 'pos'))

write.table(dat, 'data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', col.names = T, row.names = F, quote = F)
