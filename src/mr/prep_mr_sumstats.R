# Prep the summary stats for MR analysis
library(vroom)
library(dplyr)
library(purrr)
library(tidyr)

file_list = readLines('cand_compound_for_mr.txt')

# Make the colnames consistent and pull out relevant columns
clean_dat = function(dat) {
	if('P' %in% colnames(dat)) {
		colnames(dat) = c("chr", "pos", "other_allele", "effect_allele", "P", "beta", "se", "maf", "cpid")
	} else if ('p_value' %in% colnames(dat)) {
		colnames(dat) = c("snp", "P", "chr", "pos", "effect_allele", "other_allele", "maf", "beta", "se")
		dat$cpid = paste(dat$chr, dat$pos, sep = '_')
	} else if ('P-value' %in% colnames(dat)) {
		colnames(dat) = c("cpid", "effect_allele", "other_allele", "maf", "FreqSE", "MinFreq", "MaxFreq", "beta", "se", "P", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal")
		dat = dat %>% separate(cpid, into = c('chr', 'pos'), sep = '_', remove = F)
	}
	dat = dat %>% select(cpid, chr, pos, effect_allele, other_allele, maf, beta, se, P)
	return(dat)
}

# Subset the data to the variants (and proxies) that were P <= 1e-6 in any of
# the metQTL/gout GWAS for follow-up MR:
subset_dat = function(dat, cpid_list) {
	tmp = clean_dat(dat)
	res = tmp %>% filter(cpid %in% cpid_list)
	return(res)
}

# First load in all the P <= 1e-6 cpid:
cpid_list = readLines('data/ld_proxy/all_sig_variants.txt')

# Load in all of the LD-proxies:
proxy = read.delim('results/proxy/1kgp_proxy.ld', sep = '', header = T, stringsAsFactors = F)
proxy$cpid1 = paste(proxy$CHR_A, proxy$BP_A, sep = '_')
proxy$cpid2 = paste(proxy$CHR_B, proxy$BP_B, sep = '_')

proxy_list = unique(c(proxy$cpid1, proxy$cpid2))

# Combie the lists:
cpid_list = unique(c(cpid_list, proxy_list))

# Pull out the variants from all the metQTL of interest and write it out:
dat_list = map(file_list, ~ write.table(subset_dat(vroom(.x), cpid_list), gsub('.gz', '', gsub('.*/', 'data/mr/', .x)), sep = '\t', col.names = T, row.names = F, quote = F))

# Do the same for the gout GWAS:
gout = vroom('/Volumes/archive/merrimanlab/copied/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt')
gout_dat = gout %>% filter(cpid %in% cpid_list) %>% select(SNP, cpid, CHR:BP, minor:MAF, effect:P)
colnames(gout_dat) = c("SNP", "cpid", "chr", "pos", "effect_allele", "other_allele", "maf", "beta", "se", "P")

write.table(gout_dat, 'data/mr/major_gout_gwas.txt', sep = '\t', col.names = T, row.names = F, quote = F)

