library(dplyr)
library(tidyr)

dat = read.table('data/eqtl/eqtl_results.effects_aligned.final.pp0.5.txt', sep = '\t', header = T, stringsAsFactors = F)

all = dat %>% filter(locus.type == 'cis') %>% select(gtex_id:proxy.SNP, tissue:ENSG) %>% distinct
full = dat %>% filter(cohort == 'full') %>% filter(locus.type == 'cis') %>% select(gtex_id:proxy.SNP, tissue:ENSG) %>% distinct
male = dat %>% filter(cohort == 'male') %>% filter(locus.type == 'cis') %>% select(gtex_id:proxy.SNP, tissue:ENSG) %>% distinct

all_snp = unique(c(all$SNP, all$proxy.SNP))
full_snp = unique(c(full$SNP, full$proxy.SNP))
male_snp = unique(c(male$SNP, male$proxy.SNP))

all_snp = all_snp[which(!(is.na(all_snp)))]
male_snp = male_snp[which(!(is.na(male_snp)))]
full_snp = full_snp[which(!(is.na(full_snp)))]

write.table(all, 'data/eqtl/all_cis_eqtl.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(full, 'data/eqtl/full_cis_eqtl.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(male, 'data/eqtl/male_cis_eqtl.txt', sep = '\t', col.names = T, row.names = F, quote = F)

writeLines(all_snp, 'data/eqtl/all_cis_eqtl.snp.txt')
writeLines(full_snp, 'data/eqtl/full_cis_eqtl.snp.txt')
writeLines(male_snp, 'data/eqtl/male_cis_eqtl.snp.txt')

