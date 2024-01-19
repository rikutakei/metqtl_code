library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# setwd('..')

# Load in results
met = read.table('results/mr_results/met2gout.mr_res.txt', sep = '\t', header = T, stringsAsFactors = F)
gout = read.table('results/mr_results/gout2met.mr_res.txt', sep = '\t', header = T, stringsAsFactors = F)

# Load in phenotype info:
compound_info = read.delim('data/metqtl_pheno/full_pheno.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
compound_info$plasma_id = paste(compound_info$phenocode, compound_info$summaryStatistics.plasma, sep = '_')
compound_info$plasma_id = gsub('^0_|_0$', '', compound_info$plasma_id)
compound_info = compound_info %>% select(phenostring, summaryStatistics.urine, plasma_id) %>% pivot_longer(cols = c(summaryStatistics.urine, plasma_id), names_to = 'type', values_to = 'filename')
compound_info = compound_info %>% filter(filename != '0')
compound_info$type = ifelse(compound_info$type == 'plasma_id', 'plasma', 'urine')

met = left_join(met, compound_info) %>% select(phenostring, type, method:p, filename)
gout = left_join(gout, compound_info) %>% select(phenostring, type, method:p, filename)

# NOTE: there are less phenostrings compared to filenames, due to the same
# metabolite name from urine/plasma

# Plot IVW and weighted median MR results for metQTL causing gout:
met %>%
  filter(!grepl('MR', method)) %>%
  # filter(p <= 0.05 / length(unique(filename))) %>%
  filter(p <= 0.05) %>%
  mutate(phenostring = str_to_title(phenostring), type = str_to_title(type)) %>%
  ggplot(aes(x = estimate, y = phenostring, xmin = lower, xmax = upper)) +
  theme_bw() +
  ylab(label = '') +
  xlab(label = 'MR estimate') +
  facet_grid(method ~ type) +
  geom_vline(aes(xintercept = 0), linewidth = 0.25) +
  scale_y_discrete(limits=rev) +
  labs(title = 'Causality of metQTL on gout') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_pointrange(size = 0.1, linewidth = 0.25)

ggsave('results/mr_results/met2gout.mr.png', width = 12, height = 8)

# Plot MR-egger results (to check for pleiotropy):
met %>%
  filter(method == 'MR-Egger Intercept') %>%
  filter(p <= 0.05) %>%
  mutate(phenostring = str_to_title(phenostring), type = str_to_title(type)) %>%
  ggplot(aes(x = estimate, y = phenostring, xmin = lower, xmax = upper)) +
  theme_bw() +
  ylab(label = '') +
  xlab(label = 'MR estimate') +
  facet_grid(method ~ type) +
  geom_vline(aes(xintercept = 0), linewidth = 0.25) +
  scale_y_discrete(limits=rev) +
  labs(title = 'MR-Egger intercept of metQTL on gout') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_pointrange(size = 0.1, linewidth = 0.25)

ggsave('results/mr_results/met2gout.mr_egger.png', width = 12, height = 6)

# Plot IVW and weighted median MR results for gout causing metQTL:
gout %>%
  filter(!grepl('MR', method)) %>%
  filter(p <= 0.05) %>%
  mutate(phenostring = str_to_title(phenostring), type = str_to_title(type)) %>%
  ggplot(aes(x = estimate, y = phenostring, xmin = lower, xmax = upper)) +
  theme_bw() +
  ylab(label = '') +
  xlab(label = 'MR estimate') +
  facet_grid(method ~ type) +
  geom_vline(aes(xintercept = 0), linewidth = 0.25) +
  scale_y_discrete(limits=rev) +
  labs(title = 'Causality of gout on metQTL') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_pointrange(size = 0.1, linewidth = 0.25)

ggsave('results/mr_results/gout2met.mr.png', width = 12, height = 8)

# Plot MR-egger results (to check for pleiotropy):
gout %>%
  filter(method == 'MR-Egger Intercept') %>%
  filter(p <= 0.05) %>%
  mutate(phenostring = str_to_title(phenostring), type = str_to_title(type)) %>%
  ggplot(aes(x = estimate, y = phenostring, xmin = lower, xmax = upper)) +
  theme_bw() +
  ylab(label = '') +
  xlab(label = 'MR estimate') +
  facet_grid(method ~ type) +
  geom_vline(aes(xintercept = 0), linewidth = 0.25) +
  scale_y_discrete(limits=rev) +
  labs(title = 'MR-Egger intercept of gout on metQTL') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_pointrange(size = 0.1, linewidth = 0.25)

ggsave('results/mr_results/gout2met.mr_egger.png', width = 12, height = 6)

# Make a list of significant MR metabolites for metabolite-based pathway
# enrichment analysis
met$exposure = 'metQTL'
gout$exposure = 'gout'

res = rbind(met, gout)
res$method = ifelse(res$method == 'Weighted Median', 'WM', res$method)
res$method = ifelse(res$method == 'MR-Egger', 'egger', res$method)
res$method = ifelse(res$method == 'MR-Egger Intercept', 'egger_int', res$method)

res = res %>% select(phenostring, type, method, p, exposure) %>% pivot_wider(id_cols = c('phenostring', 'type', 'exposure'), names_from = 'method', names_prefix = 'p.', values_from = 'p')

write.table(res, 'results/mr_results/mr_res.p_summary.txt', sep = '\t', col.names = T, row.names = F, quote = F)
