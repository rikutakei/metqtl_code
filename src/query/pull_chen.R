# Pull out lead variant +/-500kb from the METSIM summary stats
#
# Rscript src/query/pull_metsim.R data/metsim_data/C98.tsv
library(vroom)
library(dplyr)
library(tidyr)
library(purrr)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

snp_list = read.table('data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', header = T, stringsAsFactors = F)
snp_list = unique(snp_list$SNP)

# Pull out all the variants in the region based on all the variants present in
# full/male/female data sets. In the future when full/male/female data are used
# for coloc, this data will contain the most complete set of variants.
pull_region = function(data, variant, out_name) {
	gwas_file = gsub('SNP', variant, 'results/gwas/SEX/SNP.gwas.txt')
	gwas_file = map(c('full', 'female', 'male'), ~ gsub('SEX', .x, gwas_file))
	tmp_dat = map_dfr(gwas_file, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
	var_list = unique(tmp_dat$SNP)
	res = data %>% filter(SNP %in% var_list)
	write.table(res, out_name, sep = '\t', col.names = T, row.names = F, quote = F)
}

compound = gsub('.*/', '', args[1])
compound = gsub('_buildGRCh38.tsv', '', compound)
out_name = gsub('COM', compound, 'results/metqtl/chen/COM/SNP.metqtl.txt')
out_name = map_chr(snp_list, ~ gsub('SNP', .x, out_name))

map2(snp_list, out_name, ~ pull_region(dat, .x, .y))

