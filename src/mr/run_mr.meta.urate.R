# MR analysis of the 29 metabolites that showed causal evidence for gout, and
# test whether they are also causal for urate
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

# Read in the list of 29 metabolites:
sig_mr_met = readLines('src/mr/sig_mr_metabolites_n29.txt')

# Don't need urate:
sig_mr_met = sig_mr_met[which(sig_mr_met != 'urate')]
sig_mr_met = sig_mr_met[which(sig_mr_met != 'major_urate')]

# Only keep the 28 metabolites:
sumstats = sumstats[which(names(sumstats) %in% c(sig_mr_met))]

# Load in Major urate data:
urate = read.table('data/meta_mr/major_urate_meta1.clean.mr.txt', sep = '\t', header = T, stringsAsFactors = F)

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

# Match the variants with the urate data (some variants are missing in urate,
# compared to urate data):
snp_list = map(snp_list, ~ .x[which(.x %in% urate$SNP)])

# Keep summary stats that had >2 variants (i.e. those that have a SNP-list)
ind = which(map_dbl(snp_list, length) > 2)
snp_list = snp_list[ind]

common = names(sumstats)[which(names(sumstats) %in% names(snp_list))]

snp_list = snp_list[common]
sumstats = sumstats[common]

################################################################################
# metQTL MR (metabolite -> urate)

# For each metabolite, pull out the summary stats for those variants from the
# relevant metabolite and urate, and align the effect to the metabolite data
met_dat = map2(sumstats, snp_list, ~ .x %>% filter(SNP %in% .y))

urate_dat = map(snp_list, ~ urate %>% filter(SNP %in% .x))

# Make sure the summary stats are aligned for the risk allele in the relevant
# metabolite data:

# Flip the non-risk effects in metabolite summary stats:
ind = map(met_dat, ~ which(.x$effect < 0))
met_dat = map2(met_dat, ind, ~ align_effect(.x, .y))

# Now align the urate summary stats based on the metabolite data:
ind_list = map2(met_dat, urate_dat, ~ which(.x$effect_allele != .y$effect_allele))
urate_dat = map2(urate_dat, ind_list, ~ align_effect(.x, .y))

# Check:
ind_list = map2(met_dat, urate_dat, ~ which(.x$effect_allele == .y$effect_allele))
all(map_dbl(met_dat, nrow) == map_dbl(ind_list, length))

# Check the difference in MAF between metabolite data and urate
maf_diff = map2(met_dat, urate_dat, ~ .y %>% mutate(maf_diff = .y$MAF - .x$MAF))
diff_ind = map(maf_diff, ~ which(abs(.x$maf_diff) >= 0.5))

# Flip the alleles for those data sets that have any variants with MAF
# difference >0.5 (mainly METSIM)
ind = which(map_dbl(diff_ind, length) > 0)
if (length(ind) > 0) {
  urate_dat[ind] = map2(urate_dat[ind], diff_ind[ind], ~ align_effect(.x, .y, flip = T))
}

# Check:
maf_diff = map2(met_dat, urate_dat, ~ .y %>% mutate(maf_diff = .y$MAF - .x$MAF))
all(!map_lgl(maf_diff, ~ any(which(abs(.x$maf_diff) >= 0.5))))

# The data should be ready for MR now, so run metabolite -> urate MR:
met_mr = map2(met_dat, urate_dat, ~ run_mr(.x, .y))

# Add metabolite info to the mr result before saving
met_mr = map2_dfr(met_mr, names(met_mr), ~ .x %>% mutate(exposure = .y, outcome = 'urate') %>% select(exposure, outcome, method:p))

# Save results
write.table(met_mr, 'results/mr_results/met2urate.mr_res.meta_validation.txt', sep = '\t', col.names = T, row.names = F, quote = F)

################################################################################
# Now run urate -> metabolite MR
urate_snps = readLines('data/meta_mr/major_urate.mr_snps.txt')

# Pull out urate variants that are present in the metabolite data sets
met_dat = map(sumstats, ~ .x %>% filter(SNP %in% urate_snps))
met_dat = met_dat[which(names(met_dat) != 'urate')]
urate_dat = urate %>% filter(SNP %in% met_dat[[1]]$SNP)

# One SNP is inconsistently present in metQTL data:
not_present = map(met_dat, ~ urate_dat$SNP[!(urate_dat$SNP %in% .x$SNP)])
not_present = not_present %>% unlist %>% unique

met_dat = map(met_dat, ~ .x %>% filter(!(SNP %in% not_present)))
urate_dat = urate_dat %>% filter(!(SNP %in% not_present))

# Make sure the SNP orders are all the same
all(map_lgl(met_dat, ~ all(.x$SNP == urate_dat$SNP)))

# Align to the urate risk allele:
ind = which(urate_dat$effect < 0)
urate_dat = align_effect(urate_dat, ind)

# Now align the metabolite summary stats based on the urate data:
ind_list = map(met_dat, ~ which(.x$effect_allele != urate_dat$effect_allele))
met_dat = map2(met_dat, ind_list, ~ align_effect(.x, .y))

# Check:
ind_list = map(met_dat, ~ which(.x$effect_allele == urate_dat$effect_allele))
all(map_dbl(met_dat, nrow) == map_dbl(ind_list, length))

# Check the difference in MAF between metabolite data and urate
maf_diff = map(met_dat, ~ urate_dat %>% mutate(maf_diff = urate_dat$MAF - .x$MAF))

# If any of the variants have MAF difference >= 0.5 in any of the metabolite
# data set, then "flip" the allele to the correct effect allele
if (all(map_lgl(maf_diff, ~ any(which(abs(.x$maf_diff) >= 0.5))))) {
  diff_ind = map(maf_diff, ~ which(abs(.x$maf_diff) >= 0.5))

  # Flip the alleles for those data sets that have any variants with MAF
  # difference >0.5
  ind = which(map_dbl(diff_ind, length) > 0)
  met_dat[ind] = map(met_dat[ind], diff_ind[ind], ~ align_effect(.x, .y, flip = T))
}

# Run urate -> metabolite MR
urate_mr = map(met_dat, ~ run_mr(urate_dat, .x))

# Add metabolite info to the mr result before saving
urate_mr = map2_dfr(urate_mr, names(urate_mr), ~ .x %>% mutate(exposure = 'urate', outcome = .y) %>% select(exposure, outcome, method:p))

# Save results
write.table(urate_mr, 'results/mr_results/urate2met.mr_res.meta_validation.txt', sep = '\t', col.names = T, row.names = F, quote = F)

