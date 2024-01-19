# Run coloc of GWAS vs METSIM for all cis-eQTL loci, given a metabolite code as
# input
#
# Rscript src/coloc/run_coloc.metsim.R results/metqtl/plasma/C98 full
# Rscript src/coloc/run_coloc.metsim.R results/metqtl/urine/GCST90264179 full

library(dplyr)
library(tidyr)
library(purrr)
library(coloc)

args = commandArgs(trailingOnly = TRUE)

metqtl_pheno = read.delim('data/metqtl_pheno/full_pheno.txt', header = T, stringsAsFactors = F)

if (grepl('urine', args[1])) {
	file = gsub('.*/', '', args[1])
	sample_size = metqtl_pheno$N_urine[which(metqtl_pheno$summaryStatistics.urine == file)]
} else if (grepl('plasma', args[1])) {
	file = gsub('.*/', '', args[1])
	sample_size = metqtl_pheno$N_plasma[which((metqtl_pheno$summaryStatistics.plasma == file | metqtl_pheno$phenocode == file))]
}

# Function to read in the summary stats and add unique ID based on CHR/POS and alleles
load_data = function(file_name, threshold = 0.001) {
	data = read.table(file_name, sep = '\t', header = T, stringsAsFactors = F)
	data = data %>% filter(between(MAF, threshold, (1 - threshold)))
	data_av1 = as.numeric(unlist(lapply(data$minor, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, 100))))
	data_av2 = as.numeric(unlist(lapply(data$major, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, 100))))
	data_av = data_av1 + data_av2
	data$id = paste(data$CHR, data$BP, data_av, sep = '_')
	return(data)
}

# Load in variants of interest
loci = read.table('data/gwas/indep_snps_from_paper.clean.txt', sep = '\t', header = T, stringsAsFactors = F)

met_dir = args[1]
metabolite = gsub('.*/', '', args[1])
met_files = paste(args[1], '/', unique(loci$SNP), '.metqtl.txt', sep = '')
met_list = map(met_files, ~ load_data(.x))
met_names = gsub('.*/', '', met_files)
met_names = gsub('.metqtl.txt', '', met_names)
names(met_list) = met_names

cohort = args[2]
gwas_files = paste('results/gwas/', cohort, '/', unique(loci$SNP), '.gwas.txt', sep = '')
gwas_list = map(gwas_files, ~ load_data(.x))
gwas_names = gsub('.*/', '', gwas_files)
gwas_names = gsub('.gwas.txt', '', gwas_names)
names(gwas_list) = gwas_names

# Remove those with < 100 variants
rm_met = map_lgl(met_list, ~ nrow(.x) < 100)
rm_gwas = map_lgl(gwas_list, ~ nrow(.x) < 100)

keep = names(gwas_list)[which(!(rm_met | rm_gwas))]

met_list = met_list[keep]
gwas_list = gwas_list[keep]

# Match up the files and find variants that are common
#
# NOTE: there may be duplicated variants/SNPs in metabolite data, so take care
# of it (I've also done it for GWAS, just in case)
get_common_snps = function(gwas, met) {
	# GWAS
	dup_gwas = gwas$id[duplicated(gwas$id)]
	gwas_id = gwas$id[which(!(gwas$id %in% dup_gwas))]
	# metQTL
	dup_met = met$id[duplicated(met$id)]
	met_id = met$id[which(!(met$id %in% dup_met))]
	# Common
	common = table(c(gwas_id, met_id)) == 2
	id_list = names(table(c(gwas_id, met_id)))[common]
	return(id_list)
}

common_snps = map(names(gwas_list), ~ get_common_snps(gwas_list[[.x]], met_list[[.x]]))
names(common_snps) = names(gwas_list)

# Take common SNPs
met_names = names(met_list)
met_list = map(names(met_list), ~ met_list[[.x]] %>% filter(id %in% common_snps[[.x]]))
names(met_list) = met_names

gwas_names = names(gwas_list)
gwas_list = map(names(gwas_list), ~ gwas_list[[.x]] %>% filter(id %in% common_snps[[.x]]))
names(gwas_list) = gwas_names

# If any of the locus has < 100 variants after taking the common SNPs, exclude
# them
keep = names(met_list)[map_lgl(met_list, ~ nrow(.x) >= 100)]

met_list = met_list[keep]
gwas_list = gwas_list[keep]

# Function to align the effect alleles based on the first data set (i.e. GWAS)
#
# Since only the common variants are kept so far, only need to check the effect
# allele and decide whether it needs to be flipped or not
align_allele = function(gwas, met) {
	# Just make sure the variants are in the right order:
	gwas = gwas %>% arrange(CHR, BP)
	met = met %>% arrange(CHR, BP)
	# Now check for effect allele
	ind = which(gwas$minor != met$minor)
	if (length(ind) == 0) {
		return(met)
	}
	res = met
	tmp = res$major[ind]
	res$major[ind] = res$minor[ind]
	res$minor[ind] = tmp
	res$MAF[ind] = 1 - res$MAF[ind]
	res$effect[ind] = -res$effect[ind]
	return(res)
}

met_list = map2(gwas_list, met_list, ~ align_allele(.x, .y))

# Function to pull out relevant coloc info from summary stats
# Case proportion for EUR GWAS:
#  - full = 100661 / 2106003 = 0.0478
#  - female = 22398 / 1146092 = 0.0195
#  - male = 72799 / 926287 = 0.0808
# Defaults to case proportion from full cohort
run_coloc = function(gwas_dat, met_dat, prop_case = 0.0478, n_met) {
	gwas = list(snp = gwas_dat$id,
				beta = gwas_dat$effect,
				MAF = gwas_dat$MAF,
				N = gwas_dat$N,
				varbeta = (gwas_dat$SE ^ 2),
				pvalues = gwas_dat$P,
				type = 'cc',
				s = prop_case)

	met = list(snp = met_dat$id,
			   beta = met_dat$effect,
			   MAF = met_dat$MAF,
			   N = n_met,
			   varbeta = (met_dat$SE ^ 2),
			   pvalues = met_dat$P,
			   type = 'quant')

	results = coloc.abf(gwas, met)
	return(results$summary)
}

if(cohort == 'full') {
	prop = 0.0478
} else if(cohort == 'male') {
	prop = 0.0808
} else {
	prop = 0.0195
}

coloc_res = map_dfr(names(gwas_list), ~ run_coloc(gwas_list[[.x]], met_list[[.x]], prop_case = prop, n_met = sample_size))

coloc_res$SNP = names(gwas_list)
coloc_res = coloc_res %>% select(SNP, nsnps:PP.H4.abf)

# Save result
out_name = paste(met_dir, metabolite, sep = '/')
out_name = paste(out_name, 'coloc_res', cohort, 'txt', sep = '.')
write.table(coloc_res, out_name, sep = '\t', col.names = T, row.names = F, quote = F)
