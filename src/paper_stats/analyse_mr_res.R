library(dplyr)
library(purrr)
library(tidyr)

setwd('..')

# Load in metabolite-to-gout results:
res = read.table('results/mr_results/met2gout.mr_res.txt', sep = '\t', header = T, stringsAsFactors = F)

# Load the coloc data just for converting the metabolite ID to name
file_list = c('results/coloc_res/metsim/merged_res/coloc_merged.pp0.8.txt', 'results/coloc_res/schlosser/plasma/merged_res/coloc_merged.pp0.8.plasma.txt', 'results/coloc_res/chen/merged_res/coloc_merged.pp0.8.txt')
pheno = map_dfr(file_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F, quote = '') %>% select(compound, phenostring))
pheno = pheno %>% distinct

# Add the metabolite names:
res = left_join(res, pheno, by = c('exposure' = 'compound'))

# Take a look at the MR results using Chen data:
chen = res %>% filter(dataset == 'chen')

# Note that the result includes ALL significant metabolites (including from
# Schlosser and METSIM), so pull out only the relevant metabolites for Chen:
met = readLines('results/mr_results/chen.mr_list.specific.txt')
chen = chen %>% filter(exposure %in% met)

# Pull out metabolites that had multiple testing-corrected significant IVW or
# weighted median p-vlaue
num_test = length(unique(chen$exposure))
threshold = 0.05 / num_test

sig_ivw = chen %>% filter(method == 'IVW', p <= threshold)
sig_med = chen %>% filter(method == 'Weighted Median', p <= threshold)
sig_egger = chen %>% filter(method == 'MR-Egger Intercept', p <= threshold)

sig_chen = chen %>% filter(exposure %in% c(sig_ivw$exposure, sig_med$exposure))

# Now take a look at the gout-to-metabolite results:
gout_res = read.table('results/mr_results/gout2met.mr_res.txt', sep = '\t', header = T, stringsAsFactors = F)

gout_res = left_join(gout_res, pheno, by = c('outcome' = 'compound'))

# Pull out the significant results using Chen:
gout_chen = gout_res %>% filter(dataset == 'chen', outcome %in% chen$exposure)
threshold = 0.05 / length(unique(gout_chen$outcome))

gout_ivw = gout_chen %>% filter(method == 'IVW', p <= threshold)
gout_med = gout_chen %>% filter(method == 'Weighted Median', p <= threshold)
gout_egger = gout_chen %>% filter(method == 'MR-Egger Intercept', p <= threshold)
