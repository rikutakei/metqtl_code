library(dplyr)
library(purrr)
library(tidyr)
library(vroom)
library(MendelianRandomization)

setwd('..')

source('src/mr/mr_functions.R')

# Load in all the metQTL summary statistics for MR
file_list = list.files('data/meta_mr/', pattern = '*.mr.txt', full.names = T)

sumstats = map(file_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

# Add names to the list:
name_list = gsub('.*/', '', file_list)
name_list = gsub('_meta1.*', '', name_list)
names(sumstats) = name_list

# Only keep 141 metabolites (and gout):
overlap_met = readLines("results/met_overlap/all_plasma_overlap.txt")
overlap_met = gsub('[\\(\\)]', '', overlap_met)
overlap_met = gsub(' ', '-', overlap_met)
overlap_met = gsub(',', '_', overlap_met)

sumstats = sumstats[c("gout", overlap_met)]

# Load in the selected MR variants:
snp_file = list.files('data/meta_mr/', pattern = '*.mr_snps.txt', full.names = T)

# Don't need "all_mr_snps.txt"
ind = which(grepl('all_mr', snp_file))
snp_file = snp_file[-ind]

snp_list = map(snp_file, ~ readLines(.x))

# Add names to the list:
name_list = gsub('.*/', '', snp_file)
name_list = gsub('\\..*', '', name_list)
names(snp_list) = name_list

# Keep summary stats that had >2 variants (i.e. those that have a SNP-list)
ind = which(map_dbl(snp_list, length) > 2)
snp_list = snp_list[ind]

common = names(sumstats)[which(names(sumstats) %in% names(snp_list))]

snp_list = snp_list[common]
sumstats = sumstats[common]

################################################################################
# metQTL MR (metabolite -> gout)

# For each metabolite, pull out the summary stats for those variants from the
# relevant metabolite and gout, and align the effect to the metabolite data
met_dat = map2(sumstats, snp_list, ~ .x %>% filter(SNP %in% .y))

gout_dat = read.table('data/meta_mr/gout_meta1.clean.mr.txt', sep = '\t', header = T, stringsAsFactors = F)
colnames(gout_dat) = c('SNP', 'CHR', 'POS', 'effect_allele', 'other_allele', 'MAF', 'effect', 'SE', 'P')
gout_dat = map(snp_list, ~ gout_dat %>% filter(SNP %in% .x))

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
met_mr = map2_dfr(met_mr, names(met_mr), ~ .x %>% mutate(exposure = .y, outcome = 'gout') %>% select(exposure, outcome, method:p))

# Save results
write.table(met_mr, 'results/mr_results/met2gout.mr_res.meta_validation.txt', sep = '\t', col.names = T, row.names = F, quote = F)

################################################################################
# Now run gout -> metabolite MR
gout_dat = read.table('data/meta_mr/gout_meta1.clean.mr.txt', sep = '\t', header = T, stringsAsFactors = F)
colnames(gout_dat) = c('SNP', 'CHR', 'POS', 'effect_allele', 'other_allele', 'MAF', 'effect', 'SE', 'P')
gout_snps = readLines('data/meta_mr/gout.mr_snps.txt')

# Pull out gout variants that are present in the metabolite data sets
met_dat = map(sumstats, ~ .x %>% filter(SNP %in% gout_snps))
met_dat = met_dat[which(names(met_dat) != 'gout')]
gout_dat = gout_dat %>% filter(SNP %in% met_dat[[1]]$SNP)

# Make sure the SNP orders are all the same
all(map_lgl(met_dat, ~ all(.x$SNP == gout_dat$SNP)))

# Align to the gout risk allele:
ind = which(gout_dat$effect < 0)
gout_dat = align_effect(gout_dat, ind)

# Now align the metabolite summary stats based on the gout data:
ind_list = map(met_dat, ~ which(.x$effect_allele != gout_dat$effect_allele))
met_dat = map2(met_dat, ind_list, ~ align_effect(.x, .y))

# Check:
ind_list = map(met_dat, ~ which(.x$effect_allele == gout_dat$effect_allele))
all(map_dbl(met_dat, nrow) == map_dbl(ind_list, length))

# Check the difference in MAF between metabolite data and gout
maf_diff = map(met_dat, ~ gout_dat %>% mutate(maf_diff = gout_dat$MAF - .x$MAF))

# If any of the variants have MAF difference >= 0.5 in any of the metabolite
# data set, then "flip" the allele to the correct effect allele
if (all(map_lgl(maf_diff, ~ any(which(abs(.x$maf_diff) >= 0.5))))) {
  diff_ind = map(maf_diff, ~ which(abs(.x$maf_diff) >= 0.5))

  # Flip the alleles for those data sets that have any variants with MAF
  # difference >0.5
  ind = which(map_dbl(diff_ind, length) > 0)
  met_dat[ind] = map(met_dat[ind], diff_ind[ind], ~ align_effect(.x, .y, flip = T))
}

# Run gout -> metabolite MR
gout_mr = map(met_dat, ~ run_mr(gout_dat, .x))

# Add metabolite info to the mr result before saving
gout_mr = map2_dfr(gout_mr, names(gout_mr), ~ .x %>% mutate(exposure = 'gout', outcome = .y) %>% select(exposure, outcome, method:p))

# Save results
write.table(gout_mr, 'results/mr_results/gout2met.mr_res.meta_validation.txt', sep = '\t', col.names = T, row.names = F, quote = F)

