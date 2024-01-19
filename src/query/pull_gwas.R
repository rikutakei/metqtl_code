# Pull out lead variant +/-500kb from the gout summary stats
library(vroom)
library(dplyr)
library(tidyr)
library(purrr)

full = vroom('data/gwas/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt')
male = vroom('data/gwas/EUR_meta_male1_clean_rsid.nfiltered.biallelic.txt')
female = vroom('data/gwas/EUR_meta_female1_clean_rsid.nfiltered.biallelic.txt')

snp_list = read.table('data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', header = T, stringsAsFactors = F)
snp_list = snp_list %>% distinct(SNP, chr, pos)

pull_region = function(data, chr, pos, snp, out_name) {
	start = pos - 500000
	end = pos + 500000
	res = data %>% filter(CHR == chr) %>% filter(between(BP, start, end))
	write.table(res, out_name, sep = '\t', col.names = T, row.names = F, quote = F)
}

out_name = map(snp_list$SNP, ~ gsub('SNP', .x, 'results/gwas/SEX/SNP.gwas.txt'))

out_full = map_chr(out_name, ~ gsub('SEX', 'full', .x))
out_male = map_chr(out_name, ~ gsub('SEX', 'male', .x))
out_female = map_chr(out_name, ~ gsub('SEX', 'female', .x))

pmap(list(snp_list$chr, snp_list$pos, snp_list$SNP, out_full), ~ pull_region(full, ..1, ..2, ..3, ..4))
pmap(list(snp_list$chr, snp_list$pos, snp_list$SNP, out_male), ~ pull_region(male, ..1, ..2, ..3, ..4))
pmap(list(snp_list$chr, snp_list$pos, snp_list$SNP, out_female), ~ pull_region(female, ..1, ..2, ..3, ..4))

