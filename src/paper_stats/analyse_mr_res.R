library(dplyr)
library(purrr)
library(tidyr)

setwd('..')

# Load in metabolite-to-gout results:
res = read.table('results/mr_results/met2gout.mr_res.meta_validation.txt', sep = '\t', header = T, stringsAsFactors = F)

# Remove gout:
res = res %>% filter(exposure != 'gout')

# Pull out metabolites that had multiple testing-corrected significant IVW or
# weighted median p-vlaue
num_test = length(unique(res$exposure))
threshold = 0.05 / num_test

sig_ivw = res %>% filter(method == 'IVW', p <= threshold)
sig_med = res %>% filter(method == 'Weighted Median', p <= threshold)
sig_egger = res %>% filter(method == 'MR-Egger Intercept', p <= threshold)

sig_res = res %>% filter(exposure %in% c(sig_ivw$exposure, sig_med$exposure))

length(unique(sig_res$exposure))

write.table(sig_res, 'results/mr_results/met2gout.mr_res.meta_validation.sig.txt', sep = '\t', row.names = F, col.names = T, quote = F)

# Which of the significant metabolites had pleiotropy:
sig_res_egger = sig_res %>% filter(exposure %in% c(sig_egger$exposure))

length(unique(sig_res_egger$exposure))

write.table(sig_res_egger, 'results/mr_results/met2gout.mr_res.meta_validation.sig_egger.txt', sep = '\t', row.names = F, col.names = T, quote = F)

################################################################################
# Now take a look at the gout-to-metabolite results:
gout_res = read.table('results/mr_results/gout2met.mr_res.meta_validation.txt', sep = '\t', header = T, stringsAsFactors = F)

# Pull out the significant metQTL MR results:
threshold = 0.05 / length(unique(gout_res$outcome))

gout_ivw = gout_res %>% filter(method == 'IVW', p <= threshold)
gout_med = gout_res %>% filter(method == 'Weighted Median', p <= threshold)
gout_egger = gout_res %>% filter(method == 'MR-Egger Intercept', p <= threshold)

sig_gout = gout_res %>% filter(outcome %in% c(gout_ivw$outcome, gout_med$outcome))

length(unique(sig_gout$exposure))

write.table(sig_gout, 'results/mr_results/gout2met.mr_res.meta_validation.sig.txt', sep = '\t', row.names = F, col.names = T, quote = F)

################################################################################
# Load in metabolite-to-urate results:
res = read.table('results/mr_results/met2urate.mr_res.meta_validation.txt', sep = '\t', header = T, stringsAsFactors = F)

# Remove urate:
res = res %>% filter(exposure != 'urate')

# Pull out metabolites that had multiple testing-corrected significant IVW or
# weighted median p-vlaue
num_test = length(unique(res$exposure))
threshold = 0.05 / num_test

sig_ivw = res %>% filter(method == 'IVW', p <= threshold)
sig_med = res %>% filter(method == 'Weighted Median', p <= threshold)
sig_egger = res %>% filter(method == 'MR-Egger Intercept', p <= threshold)

sig_res = res %>% filter(exposure %in% c(sig_ivw$exposure, sig_med$exposure))

write.table(sig_res, 'results/mr_results/met2urate.mr_res.meta_validation.sig.txt', sep = '\t', row.names = F, col.names = T, quote = F)

# Which of the significant metabolites had pleiotropy:
sig_res_egger = sig_res %>% filter(exposure %in% c(sig_egger$exposure))

length(unique(sig_res_egger$exposure))

write.table(sig_res_egger, 'results/mr_results/met2urate.mr_res.meta_validation.sig_egger.txt', sep = '\t', row.names = F, col.names = T, quote = F)

# What were the non-significant metabolites:
non_sig = res %>% filter(!(exposure %in% sig_res$exposure))

unique(non_sig$exposure)

################################################################################
# urate-to-metabolite results:
urate_res = read.table('results/mr_results/urate2met.mr_res.meta_validation.txt', sep = '\t', header = T, stringsAsFactors = F)

# Pull out the significant metQTL MR results:
threshold = 0.05 / length(unique(urate_res$outcome))

urate_ivw = urate_res %>% filter(method == 'IVW', p <= threshold)
urate_med = urate_res %>% filter(method == 'Weighted Median', p <= threshold)
urate_egger = urate_res %>% filter(method == 'MR-Egger Intercept', p <= threshold)

sig_urate = urate_res %>% filter(outcome %in% c(urate_ivw$outcome, urate_med$outcome))

length(unique(sig_urate$exposure))

write.table(sig_urate, 'results/mr_results/urate2met.mr_res.meta_validation.sig.txt', sep = '\t', row.names = F, col.names = T, quote = F)

citation()

