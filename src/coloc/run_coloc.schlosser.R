# Run coloc of GWAS vs METSIM for all cis-eQTL loci, given a metabolite code as
# input
#
# Rscript src/coloc/run_coloc.schlosser.R GCST90264176 full <plasma/urine>

library(dplyr)
library(tidyr)
library(purrr)
library(coloc)

args = commandArgs(trailingOnly = TRUE)

type = args[3]

metabolite = args[1]
met_dir = paste('results/metqtl/schlosser/', type, sep = '')
met_dir = paste(met_dir, metabolite, sep = '/')
met_files = list.files(met_dir, pattern = 'metqtl.txt', full.names = T)
met_list = map(met_files, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
met_names = gsub('.*/', '', met_files)
met_names = gsub('.metqtl.txt', '', met_names)
names(met_list) = met_names

# Use full data set for now:
cohort = args[2]
gwas_files = list.files(paste('results/gwas/', cohort, '/', sep = ''), pattern = 'gwas.txt', full.names = T)
gwas_list = map(gwas_files, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
gwas_names = gsub('.*/', '', gwas_files)
gwas_names = gsub('.gwas.txt', '', gwas_names)
names(gwas_list) = gwas_names

# Match up the files and find variants that are common
# NOTE: there may be duplicated variants/SNPs in metabolite data, so take care
# of it (I've also done it for GWAS, just in case)
get_common_snps = function(gwas, met) {
	gwas_snp = gwas$SNP
	dup_gwas = gwas_snp[duplicated(gwas_snp)]
	gwas_snp = gwas_snp[which(!(gwas_snp %in% dup_gwas))]
	met_snp = met$SNP
	dup_met = met_snp[duplicated(met_snp)]
	met_snp = met_snp[which(!(met_snp %in% dup_met))]
	common = table(c(gwas_snp, met_snp)) == 2
	snp_list = names(table(c(gwas_snp, met_snp)))[common]
	return(snp_list)
}

common_snps = map(names(gwas_list), ~ get_common_snps(gwas_list[[.x]], met_list[[.x]]))
names(common_snps) = names(gwas_list)

met_list = map(names(common_snps), ~ met_list[[.x]] %>% filter(SNP %in% common_snps[[.x]]))
gwas_list = map(names(common_snps), ~ gwas_list[[.x]] %>% filter(SNP %in% common_snps[[.x]]))

names(met_list) = met_names
names(gwas_list) = gwas_names

# If any of the locus has < 100 variants, exclude them
keep = map_lgl(met_list, ~ nrow(.x) >= 100)

met_list = met_list[keep]
gwas_list = gwas_list[keep]

# Function to pull out relevant coloc info from summary stats
# Case proportion for EUR GWAS:
#  - full = 100661 / 2106003 = 0.0478
#  - female = 22398 / 1146092 = 0.0195
#  - male = 72799 / 926287 = 0.0808
# Defaults to case proportion from full cohort
run_coloc = function(gwas_dat, met_dat, prop_case = 0.0478) {
	gwas = list(snp = gwas_dat$SNP,
				beta = gwas_dat$effect,
				MAF = gwas_dat$MAF,
				N = gwas_dat$N,
				varbeta = (gwas_dat$SE ^ 2),
				pvalues = gwas_dat$P,
				type = 'cc',
				s = prop_case)

	met = list(snp = met_dat$SNP,
			   beta = met_dat$effect,
			   MAF = met_dat$MAF,
			   N = 5023,
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

coloc_res = map_dfr(names(gwas_list), ~ run_coloc(gwas_list[[.x]], met_list[[.x]], prop_case = prop))

coloc_res$SNP = names(gwas_list)
coloc_res = coloc_res %>% select(SNP, nsnps:PP.H4.abf)

# Save result
out_name = paste('results/coloc_res/schlosser', type, cohort, metabolite, sep = '/')
out_name = paste(out_name, 'coloc_res', cohort, 'txt', sep = '.')
write.table(coloc_res, out_name, sep = '\t', col.names = T, row.names = F, quote = F)
