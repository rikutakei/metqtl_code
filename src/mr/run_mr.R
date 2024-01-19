library(dplyr)
library(purrr)
library(tidyr)
library(vroom)
library(MendelianRandomization)

# setwd('..')

# Load in all the metQTL summary statistics for MR
file_list = list.files('data/mr/sum_stats/', full.names = T)
file_list = file_list[which(!grepl('major', file_list))]

sumstats = map(file_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
sumstats = map(sumstats, ~ .x %>% mutate(chr = as.numeric(chr), pos = as.numeric(pos), effect_allele = toupper(effect_allele), other_allele = toupper(other_allele)) %>% filter(!is.na(pos)))

# Load gout summary stats separately
gout = read.table('data/mr/sum_stats/major_gout_gwas.txt', sep = '\t', header = T, stringsAsFactors = F)
gout_sig = gout %>% filter(P <= 1e-6)

# Load in some additional data (lead variants and proxies)
all_lead = readLines('data/mr/all_sig_variants.txt')

proxy = read.table('data/mr/1kgp_proxy.ld', sep = '', header = T, stringsAsFactors = F)
proxy = proxy %>% filter(SNP_A != SNP_B) %>% mutate(cpid1 = paste(CHR_A, BP_A, sep = '_'), cpid2 = paste(CHR_B, BP_B, sep = '_'))

# Little bit complicated, but I need to figure out what the lead variants for
# each metQTL data are. It's strightforward if the data only had the lead
# variants for that metQTL, but it contains ALL lead variants from all metQTL,
# plus the LD-proxies.
#
# The approach I'm taking is to pull out all the lead variants first, then
# filter for P <= 1e-6. These variants will be the lead variants for the metQTL
# data. At a later stage when I merge the gout data with metQTL data, I will
# pull out proxy data from the data.
lead_list = map(sumstats, ~ .x %>% filter(cpid %in% all_lead, P <= 1e-6))

# There are some that only have a few variants with P <= 1e-6 - these are
# likely from plasma/urine counterparts where it was colocalised in plasma but
# not urine (I pulled out the summary stats for both plasma and urine, if it
# was available)
map_dbl(lead_list, nrow)

# Assign a locus to significant SNPs by checking how close they are to one
# another.
# If the variants are <=50kb away, I'll consider it as the same locus.
assign_locus = function(dat, dist = 50000) {
  res = dat %>% arrange(chr, pos)
  res$locus = 0
  locus = 1
  res$locus[1] = locus
  if (nrow(res) == 1) {
    return(res)
  }
  for (i in 2:nrow(res)) {
    if (res$chr[i] == res$chr[i -1] & (res$pos[i] - res$pos[i - 1]) <= dist) {
      res$locus[i] = locus
      next
    }
    locus = locus + 1
    res$locus[i] = locus
  }
  return(res)
}

met_sig = map(lead_list, ~ assign_locus(.x))
gout_sig = assign_locus(gout_sig)

# NOTE: gout_sig and met_sig will be used for MR for gout -> metQTL and
# metQTL -> gout, respectively.
# i.e. met_sig will be used to test for causal relationship of metQTL on gout,
# and vice versa

# First do the metQTL MR

# Function to merge the GWAS and metQTL data for common variants
merge_dat = function(dat1, dat2, suf = c('.metqtl', '.gout')) {
  res = left_join(dat1, dat2, by = c('cpid', 'chr', 'pos'), suffix = suf)
  return(res)
}

comb_met = map(met_sig, ~ merge_dat(.x, gout, suf = c('.metqtl', '.gout')))
comb_gout = map(sumstats, ~ merge_dat(gout_sig, .x, suf = c('.gout', '.metqtl')))

# Function to take care of loci with no common variants by looking at the
# proxies:
check_proxy = function(merged, dat1, dat2, proxy, suf = c('.metqtl', '.gout')) {
  res = merged
  res$R2 = NA
  miss = merged %>% group_by(locus) %>% filter(all(is.na(effect_allele.metqtl)) | all(is.na(effect_allele.gout))) %>% ungroup %>% pull(locus) %>% unique
  if (length(miss) > 0) {
    for (i in miss) {
      snps = merged$cpid[which(merged$locus %in% i)]
      proxy_sub = proxy %>% filter(cpid1 %in% snps)
      if (nrow(proxy_sub) == 0) { next }
      # Pull out all the proxies from both dat1 and dat2
      proxy1 = dat1 %>% filter(cpid %in% proxy_sub$cpid2)
      proxy2 = dat2 %>% filter(cpid %in% proxy_sub$cpid2)
      if (nrow(proxy1) == 0 | nrow(proxy2) == 0) { next }
      proxy1$locus = i
      proxy_comb = merge_dat(proxy1, proxy2, suf)
      proxy_comb = proxy_sub %>% select(cpid2, R2) %>% left_join(proxy_comb, ., by = c('cpid' = 'cpid2'))
      res = rbind(res, proxy_comb)
    }
  }
  return(res)
}

comb_met = map2(comb_met, sumstats, ~ check_proxy(.x, .y, gout, proxy, suf = c('.metqtl', '.gout')))
comb_gout = map2(comb_gout, sumstats, ~ check_proxy(.x, gout, .y, proxy, suf = c('.gout', '.metqtl')))

# Now remove those without gout/metQTL data:
comb_met_clean = map(comb_met, ~ .x %>% filter(!is.na(effect_allele.gout)))
comb_gout_clean = map(comb_gout, ~ .x %>% filter(!is.na(effect_allele.metqtl)))

# Pick out appropriate lead variants for each locus
get_lead = function(dat, p_from = 'met') {
  # Pull out loci lead variants from lead variants and proxy variants (if none
  # of the lead variants at the locus was present in both data set)
  #
  # NOTE: loci with proxies are, by default, loci that didn't have variants in
  # both metQTL/gout, since that is why the proxies were added in. Therefore,
  # proxy and non-proxy list should not have overlapping loci
  non_proxy = dat %>% filter(is.na(R2))
  proxy_dat = dat %>% filter(!is.na(R2))
  if (p_from == 'met') {
    non_proxy = non_proxy %>% group_by(locus) %>% arrange(P.metqtl) %>% slice(1) %>% ungroup
    proxy_dat = proxy_dat %>% group_by(locus) %>% arrange(P.metqtl) %>% slice(1) %>% ungroup
  } else {
    non_proxy = non_proxy %>% group_by(locus) %>% arrange(P.gout) %>% slice(1) %>% ungroup
    proxy_dat = proxy_dat %>% group_by(locus) %>% arrange(P.gout) %>% slice(1) %>% ungroup
  }
  res = rbind(non_proxy, proxy_dat) %>% arrange(locus)
  return(res)
}

mr_met = map(comb_met_clean, ~ get_lead(.x, p_from = 'met'))
mr_gout = map(comb_gout_clean, ~ get_lead(.x, p_from = 'gout'))

# Function to check for LD between the chosen lead variants
check_ld = function(dat, proxy, p_from = 'met') {
  if (p_from == 'met') {
    dat = dat %>% arrange(P.metqtl)
  } else {
    dat = dat %>% arrange(P.gout)
  }
  proxy_sub = proxy %>% filter(cpid1 %in% dat$cpid)
  res = data.frame()
  for (i in 1:nrow(dat)) {
    ld_var = proxy_sub %>% filter(cpid1 == dat$cpid[i]) %>% pull(cpid2)
    ld_var = c(dat$cpid[i], ld_var)
    ld_dat = dat %>% filter(cpid %in% ld_var) %>% slice(1)
    res = rbind(res, ld_dat)
  }
  res = res %>% distinct %>% arrange(chr, pos) %>% as.data.frame
  return(res)
}

mr_met = map(mr_met, ~ check_ld(.x, proxy, p_from = 'met'))
mr_gout = map(mr_gout, ~ check_ld(.x, proxy, p_from = 'gout'))

# The data should be ready for MR, so now it's just a matter of adjusting the
# data a little bit for MR
#
# Align the alleles so the effects are for the same allele:
# NOTE: MAF is ignored, since MR doesn't require it
align_effect = function(dat) {
  ind = which(dat$effect_allele.metqtl != dat$effect_allele.gout)
  if (length(ind) == 0) { return(dat) }
  res = dat
  tmp = res$effect_allele.metqtl[ind]
  res$effect_allele.metqtl[ind] = res$other_allele.metqtl[ind]
  res$other_allele.metqtl[ind] = tmp
  res$beta.metqtl[ind] = -res$beta.metqtl[ind]
  return(res)
}

mr_met_clean = map(mr_met, ~ align_effect(.x))
mr_gout_clean = map(mr_gout, ~ align_effect(.x))

# Run Mendelian randomisation (X = exposure, Y = outcome)
run_mr = function(dat, exposure = 'met') {
  if (exposure == 'met') {
    input = mr_input(bx = dat$beta.metqtl, bxse = dat$se.metqtl, by = dat$beta.gout, byse = dat$se.gout)
  } else {
    input = mr_input(bx = dat$beta.gout, bxse = dat$se.gout, by = dat$beta.metqtl, byse = dat$se.metqtl)
  }
  if (length(input@snps) <= 2) { return("not enough variants")}
  # Run various MR methods:
  ivw = mr_ivw(input)
  med = mr_median(input, weighting = 'weighted', iterations = 10000, distribution = 'normal')
  egger = mr_egger(input)
  # Pull out relevant stats:
  ivw_res = data.frame(method = 'IVW',
                       estimate = ivw@Estimate,
                       se = ivw@StdError,
                       lower = ivw@CILower,
                       upper = ivw@CIUpper,
                       p = ivw@Pvalue)
  med_res = data.frame(method = 'Weighted Median',
                       estimate = med@Estimate,
                       se = med@StdError,
                       lower = med@CILower,
                       upper = med@CIUpper,
                       p = med@Pvalue)
  egger_res = data.frame(method = c('MR-Egger', 'MR-Egger Intercept'),
                         estimate = c(egger@Estimate, egger@Intercept),
                         se = c(egger@StdError.Est, egger@StdError.Int),
                         lower = c(egger@CILower.Est, egger@CILower.Int),
                         upper = c(egger@CIUpper.Est, egger@CIUpper.Int),
                         p = c(egger@Pvalue.Est, egger@Pvalue.Int))
  res = rbind(ivw_res, med_res, egger_res)
  return(res)
}

met_mr_res = map(mr_met_clean, ~ run_mr(.x, exposure = 'met'))
gout_mr_res = map(mr_gout_clean, ~ run_mr(.x, exposure = 'gout'))

# Assign original filenames:
name_list = gsub('.*/', '', file_list)
name_list = gsub('\\..*', '', name_list)
name_list = gsub('_build.*', '', name_list)

names(met_mr_res) = name_list
names(gout_mr_res) = name_list

# Keep only those that have MR results:
ind = which(map_lgl(met_mr_res, ~ is.data.frame(.x)))
met_mr_res = met_mr_res[ind]

ind = which(map_lgl(gout_mr_res, ~ is.data.frame(.x)))
gout_mr_res = gout_mr_res[ind]

# Add the filename/metQTL name to the data frames
met_mr_res = map2_dfr(met_mr_res, names(met_mr_res), ~ .x %>% mutate(filename = .y))
gout_mr_res = map2_dfr(gout_mr_res, names(gout_mr_res), ~ .x %>% mutate(filename = .y))

met_mr_res = met_mr_res %>% select(filename, method:p)
gout_mr_res = gout_mr_res %>% select(filename, method:p)

# Save results
write.table(met_mr_res, 'results/mr_results/met2gout.mr_res.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(gout_mr_res, 'results/mr_results/gout2met.mr_res.txt', sep = '\t', col.names = T, row.names = F, quote = F)

