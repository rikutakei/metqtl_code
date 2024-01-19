# Pull out lead variant +/-500kb from the urine metQTL summary stats from
# Schlosser
#
# Rscript src/query/pull_urine_metqtl.R data/schlosser_metqtl/GCST90266422_buildGRCh37.tsv.gz
library(vroom)
library(dplyr)
library(tidyr)
library(purrr)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])
colnames(dat) = c('SNP', 'P', 'CHR', 'BP', 'minor', 'major', 'MAF', 'effect', 'SE')

snp_list = read.table('data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', header = T, stringsAsFactors = F)
snp_list = snp_list %>% distinct(SNP, chr, pos)

pull_region = function(data, chr, pos) {
	start = pos - 500000
	end = pos + 500000
	res = data %>% filter(CHR == chr) %>% filter(between(BP, start, end))
	return(res)
}

dat_list = map2(snp_list$chr, snp_list$pos, ~ pull_region(dat, .x, .y))

# Save regions
compound = gsub('.*/', '', args[1])
compound = gsub('_buildGRCh37.tsv.gz', '', compound)
out_name = gsub('COM', compound, 'results/metqtl/schlosser/urine/COM/SNP.metqtl.txt')
out_name = map(snp_list$SNP, ~ gsub('SNP', .x, out_name))

map2(dat_list, out_name, ~ write.table(.x, .y, sep = '\t', col.names = T, row.names = F, quote = F))

