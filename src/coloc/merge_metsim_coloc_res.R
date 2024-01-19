# Script to merge all the PPH4 from all the METSIM coloc results
library(dplyr)
library(tidyr)
library(purrr)

full_list = list.files('results/coloc_res/metsim/full', full.names = T)
compounds = gsub('.*/', '', full_list)
compounds = gsub('\\..*', '', compounds)
male_list = gsub('full', 'male', full_list)
female_list = gsub('full', 'female', full_list)

res_full = map2_dfr(full_list, compounds, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(compound = .y, cohort = 'full'))
res_male = map2_dfr(male_list, compounds, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(compound = .y, cohort = 'male'))
res_female = map2_dfr(female_list, compounds, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(compound = .y, cohort = 'female'))

res_list = map_dfr(list(res_full, res_male, res_female), ~ .x)

# Add loci information - currently looking only at EUR loci/SNPs
loci = read.table('data/gwas//COMBINED_loci_summary_updated_9Dec2022_withLD_withBroad.txt', sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% filter(EUR.full != '') %>% select(locus.full, EUR.full, EUR.male, EUR.female)
loci = loci %>% pivot_longer(cols = c('EUR.full', 'EUR.male', 'EUR.female'), names_to = 'cohort', values_to = 'SNP')
loci = loci %>% mutate(cohort = gsub('EUR.', '', cohort)) %>% filter(SNP != '')
colnames(loci)[1] = 'locus'
loci$locus = gsub('_[[:alpha:]]', '', loci$locus)

res_list = res_list %>% left_join(., loci, by = 'SNP', suffix = c('.coloc', '.locus'), relationship = 'many-to-many') %>% select(locus, contains('cohort'), compound, SNP:PP.H4.abf)

# Add compound name:
compound_info = read.table('data/metsim_data/phenotypes.tsv', sep = '\t', header = T, stringsAsFactors = F, quote = '')
compound_info = compound_info %>% select(phenocode:phenostring)

res_list = res_list %>% left_join(., compound_info, by = c('compound' = 'phenocode')) %>% select(locus:compound, category:phenostring, SNP:PP.H4.abf)

# Find SNP-compound pair that had significant PPH4 (at H4 >= 0.8)
res_sig = res_list %>% filter(PP.H4.abf >= 0.8)

# For now, focus on results that have matching cohort (i.e. coloc that used the
# full cohort GWAS summary stats will only focus on SNPs from full cohort)
res_sig = res_sig %>% filter(cohort.coloc == cohort.locus)

# How many metabolites showed colocalisation with gout?
res_sig %>% distinct(cohort.coloc, compound) %>% pull(cohort.coloc) %>% table

# How many variants colocalised with each compound?
compound_count_var = res_sig %>% group_by(cohort.coloc, phenostring) %>% summarize(count = n(), coloc_variants = paste(SNP, collapse = '; ')) %>% arrange(desc(count)) %>% ungroup
table(compound_count_var$cohort.coloc, compound_count_var$count)

# How many loci colocalised with each compound?
compound_count_loci = res_sig %>% distinct(phenostring, locus, cohort.coloc) %>% group_by(cohort.coloc, phenostring) %>% summarize(count = n(), locus = paste(locus, collapse = '; ')) %>% arrange(cohort.coloc, desc(count)) %>% ungroup
table(compound_count_loci$cohort.coloc, compound_count_loci$count)

# How many compounds colocalise at a locus?
loci_count_met = res_sig %>% distinct(phenostring, locus, cohort.coloc) %>% group_by(locus, cohort.coloc) %>% summarize(count = n(), compound_list = paste(phenostring, collapse = '; ')) %>% arrange(cohort.coloc, desc(count)) %>% ungroup
table(loci_count_met$cohort.coloc, loci_count_met$count)

map_dfr(c('full', 'male', 'female'), ~ loci_count_met %>% filter(cohort.coloc == .x) %>% pull(count) %>% summary)

# Combine the per-compound summary tables
per_comp = left_join(compound_count_loci, compound_count_var, by = c('cohort.coloc', 'phenostring'), suffix = c('.locus', '.snps')) %>% arrange(phenostring, cohort.coloc, locus, count.locus, coloc_variants, count.snps)
per_comp = left_join(per_comp, compound_info, by = c('phenostring'))
per_comp = per_comp %>% select(phenocode, phenostring, category, cohort.coloc, count.locus:coloc_variants)

# Save a bunch of files
write.table(res_list, 'results/coloc_res/metsim/merged_res/coloc_merged.txt', sep = '\t', col.names = T, row.names = F, quote = F) # NOTE: this file isn't filtered for matching coloc variant and locus cohort
write.table(res_sig, 'results/coloc_res/metsim/merged_res/coloc_merged.pp0.8.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(per_comp, 'results/coloc_res/metsim/merged_res/per_compound_summary.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(loci_count_met, 'results/coloc_res/metsim/merged_res/per_locus_summary.txt', sep = '\t', col.names = T, row.names = F, quote = F)

