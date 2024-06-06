library(dplyr)
library(purrr)
library(tidyr)
library(vroom)
library(MendelianRandomization)

# setwd('..')

source('src/mr/mr_functions.R')

# Load in all the metQTL summary statistics for MR
file_list = list.files('data/mr_dat/', pattern = '*.mr.txt', full.names = T)

sumstats = map(file_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

# Add names to the list:
name_list = gsub('.*/', '', file_list)
name_list = gsub('\\..*', '', name_list)
names(sumstats) = name_list

# Load in the selected MR variants:
snp_file = list.files('data/mr_dat/', pattern = '*.mr_snps.txt', full.names = T)

# Don't need "all_mr_snps.txt"
snp_file = snp_file[2:length(snp_file)]

snp_list = map(snp_file, ~ readLines(.x))

# Add names to the list:
name_list = gsub('.*/', '', snp_file)
name_list = gsub('\\..*', '', name_list)
names(snp_list) = name_list

################################################################################
# Gout MR (gout -> metabolite)

# Pull out gout lead variants from all the summary stats:
gout_snps = snp_list[[which(names(snp_list) == 'gout')]]

sub_list = map(sumstats, ~ .x %>% filter(SNP %in% gout_snps))

# Pull out gout summary stats and a list of summary stats without gout:
gout = sub_list[[which(names(sub_list) == 'gout')]]
other_dat = sub_list[which(names(sub_list) != 'gout')]

# Make sure the gout summary stats are aligned for the risk allele:

# Flip the non-risk effects in gout summary stats:
ind = which(gout$effect < 0)
gout = align_effect(gout, ind)

# Now align the other summary stats based on the gout data:
ind_list = map(other_dat, ~ which(.x$effect_allele != gout$effect_allele))
other_dat = map2(other_dat, ind_list, ~ align_effect(.x, .y))

# Check:
ind_list = map(other_dat, ~ which(.x$effect_allele == gout$effect_allele))
all(map_dbl(ind_list, length) == nrow(gout))

# Some data may have the effects/MAF for the wrong allele, so double-check the
# MAF with the gout data to make sure the effects are for the correct allele
#
# Check the difference in MAF between gout and the other data set - variants
# with MAF difference greater than 0.5 most likely mean that the whole data set
# is referencing the wrong allele as the effect allele
maf_diff = map(other_dat, ~ .x %>% mutate(maf_diff = .x$MAF - gout$MAF))
diff_ind = map(maf_diff, ~ which(abs(.x$maf_diff) >= 0.5))

# Flip the alleles for those data sets that have any variants with MAF
# difference >0.5 (mainly METSIM)
ind = which(map_dbl(diff_ind, length) > 0)
other_dat[ind] = map2(other_dat[ind], diff_ind[ind], ~ align_effect(.x, .y, flip = T))

# Check:
maf_diff = map(other_dat, ~ .x %>% mutate(maf_diff = .x$MAF - gout$MAF))
all(!map_lgl(maf_diff, ~ any(which(abs(.x$maf_diff) >= 0.5))))

# The data should be ready for MR now, so run gout -> metabolite MR:
gout_mr = map(other_dat, ~ run_mr(gout, .x))

# Add metabolite info to the mr result before saving
gout_mr = map2_dfr(gout_mr, names(gout_mr), ~ .x %>% mutate(exposure = 'gout', metabolite = .y) %>% select(exposure, metabolite, method:p) %>% separate(metabolite, into = c('dataset', 'outcome'), sep = '_', remove = T))

# Save results
write.table(gout_mr, 'results/mr_results/gout2met.mr_res.txt', sep = '\t', col.names = T, row.names = F, quote = F)

################################################################################
# metQTL MR (metabolite -> gout)

# Make a SNP list without gout
met_snps = snp_list[which(names(snp_list) != 'gout')]

# Remove any that has <2 SNPs
ind = which(map_dbl(met_snps, length) > 2)
met_snps = met_snps[ind]

# For each metabolite, pull out the summary stats for those variants from the
# relevant metabolite and gout, and align the effect to the metabolite data
met_dat = sumstats[names(met_snps)]
met_dat = map2(met_dat, met_snps, ~ .x %>% filter(SNP %in% .y))

gout_dat = sumstats[[which(names(sumstats) == 'gout')]]
gout_dat = map(met_snps, ~ gout_dat %>% filter(SNP %in% .x))

# Make sure the summary stats are aligned for the risk allele in the relevant
# metabolite data:

# Flip the non-risk effects in metabolite summary stats:
ind = map(met_dat, ~ which(.x$effect < 0))
met_dat = map2(met_dat, ind, ~ align_effect(.x, .y))

# Now align the gout summary stats based on the metabolite data:
ind_list = map2(met_dat, gout_dat, ~ which(.x$effect_allele != .y$effect_allele))
gout_dat = map2(gout_dat, ind_list, ~ align_effect(.x, .y))

# Check:
ind_list = map2(met_dat, gout_dat, ~ which(.x$effect_allele == .y$effect_allele))
all(map_dbl(met_dat, nrow) == map_dbl(ind_list, length))

# Check the difference in MAF between metabolite data and gout
maf_diff = map2(met_dat, gout_dat, ~ .y %>% mutate(maf_diff = .y$MAF - .x$MAF))
diff_ind = map(maf_diff, ~ which(abs(.x$maf_diff) >= 0.5))

# Flip the alleles for those data sets that have any variants with MAF
# difference >0.5 (mainly METSIM)
ind = which(map_dbl(diff_ind, length) > 0)
gout_dat[ind] = map2(gout_dat[ind], diff_ind[ind], ~ align_effect(.x, .y, flip = T))

# Check:
maf_diff = map2(met_dat, gout_dat, ~ .y %>% mutate(maf_diff = .y$MAF - .x$MAF))
all(!map_lgl(maf_diff, ~ any(which(abs(.x$maf_diff) >= 0.5))))

# The data should be ready for MR now, so run metabolite -> gout MR:
met_mr = map2(met_dat, gout_dat, ~ run_mr(.x, .y))

# Add metabolite info to the mr result before saving
met_mr = map2_dfr(met_mr, names(met_mr), ~ .x %>% mutate(exposure = .y, outcome = 'gout') %>% select(exposure, outcome, method:p) %>% separate(exposure, into = c('dataset', 'exposure'), sep = '_', remove = T))

# Save results
write.table(met_mr, 'results/mr_results/met2gout.mr_res.txt', sep = '\t', col.names = T, row.names = F, quote = F)
