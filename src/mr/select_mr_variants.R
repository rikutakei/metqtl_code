# Select variants for MR analysis based on the clump results and the LD
# information for the clump lead variants

library(dplyr)
library(tidyr)
library(purrr)

# Load the clumping results and order them by P-value (most to least significant):
clump_list = list.files('data/mr_dat/', pattern = '*.clumped', full.names = T)
clump_dat = map(clump_list, ~ read.delim(.x, sep = '', header = T, stringsAsFactors = F) %>% select(CHR:P) %>% arrange(P))
clump_name = gsub('.sumstats.*', '', clump_list)
clump_name = gsub('.*/', '', clump_name)
names(clump_dat) = clump_name

# Load the LD information:
ld_list = list.files('data/mr_dat/', pattern = '*.ld', full.names = T)
ld_dat = map(ld_list, ~ read.delim(.x, sep = '', header = T, stringsAsFactors = F))
ld_name = gsub('.sumstats.*', '', ld_list)
ld_name = gsub('.*/', '', ld_name)
names(ld_dat) = ld_name

# Function to pull out the most significant, independent variants for MR
# analysis
indep_var = function(dat, ld) {
	res = dat$SNP

	if (nrow(ld) == 0) {
		return(res)
	}

	# For each clump lead variant, remove any other variants that are in LD
	# from the MR variant list
	for (i in 1:nrow(dat)) {
		snp = dat$SNP[i]
		ld_check = ld %>% filter(SNP_A == snp)
		if (nrow(ld_check) > 0) {
			# Remove variants that are in LD at R2 > 0.01 from the list of
			# variants
			res = res[which(!(res %in% ld_check$SNP_B))]
			# Removed variants should also be removed from the LD data so that
			# the "good" variants don't get discarded later on
			ld = ld %>% filter(SNP_A %in% res)
		}
	}

	return(res)
}

# Since MR doesn't work with less than two SNPs, remove those metabolites
ind = map_dbl(clump_dat, nrow)
ind = which(ind <= 2)

name_list = names(clump_dat)[ind]

clump_dat = clump_dat[which(!(names(clump_dat) %in% name_list))]
ld_dat = ld_dat[which(!(names(ld_dat) %in% name_list))]

# Now generate the SNP list for the metabolites that meet the criteria
snp_list = map2(clump_dat, ld_dat, ~ indep_var(.x, .y))

# Save the SNP lists
out_list = paste('data/mr_dat/', names(snp_list), sep = '')
out_list = paste(out_list, '.mr_snps.txt', sep = '')

map2(snp_list, out_list, ~ writeLines(.x, .y))
