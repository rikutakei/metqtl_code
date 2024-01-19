# Script to merge all the PPH4 from all the coloc results for metQTL
library(dplyr)
library(tidyr)
library(purrr)

compound_info = read.delim('data/metqtl_pheno/full_pheno.txt', sep = '\t', header = T, stringsAsFactors = F)
compound_info = compound_info %>% select(phenostring:summaryStatistics.urine) %>% pivot_longer(cols = c(phenocode:summaryStatistics.urine), names_to = 'dummy', values_to = 'id') %>% filter(id != '0') %>% select(phenostring, id) %>% distinct

plasma_dir = list.files('results/metqtl/plasma/', full.names = T)
urine_dir = list.files('results/metqtl/urine/', full.names = T)

plasma_full = map_chr(plasma_dir, ~ list.files(.x, 'coloc_res.full.txt', full.names = T))
urine_full = map_chr(urine_dir, ~ list.files(.x, 'coloc_res.full.txt', full.names = T))

plasma_compound = gsub('.*/', '', plasma_dir)
urine_compound = gsub('.*/', '', urine_dir)

plasma_res = map2_dfr(plasma_full, plasma_compound, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% select(SNP, nsnps, PP.H4.abf) %>% mutate(compound = .y, type = 'plasma'))
urine_res = map2_dfr(urine_full, urine_compound, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% select(SNP, nsnps, PP.H4.abf) %>% mutate(compound = .y, type = 'urine'))
full_res = rbind(plasma_res, urine_res)

# Find SNP-compound pair that had significant PPH4 (at H4 >= 0.8)
full_res = full_res %>% filter(PP.H4.abf >= 0.8)

# Add loci information - currently looking only at EUR loci/SNPs
loci = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/COMBINED_loci_summary_updated_9Dec2022_withLD_withBroad.txt', sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% filter(EUR.full != '') %>% select(locus.full, EUR.full)

full_res = left_join(full_res, loci, by = c('SNP' = 'EUR.full')) %>% filter(!is.na(locus.full)) %>% select(locus.full, SNP:type)
colnames(full_res) = gsub('.full', '', colnames(full_res))

all_sig = full_res %>% arrange(SNP, compound, PP.H4.abf)
all_sig = left_join(all_sig, compound_info, by = c('compound' = 'id'))

# How many metabolites showed colocalisation with gout?
all_sig %>% distinct(type, compound) %>% pull(type) %>% table

# How many variants colocalised with each compound?
compound_count_var = all_sig %>% group_by(type, compound) %>% summarize(count = n(), coloc_variants = paste(SNP, collapse = '; ')) %>% left_join(., compound_info, by = c('compound' = 'id')) %>% arrange(type, desc(count)) %>% ungroup

# How many loci colocalised with each compound?
compound_count_loci = all_sig %>% mutate(locus_tmp = gsub('_[A-Z]', '', locus)) %>% distinct(compound, locus_tmp, type) %>% group_by(type, compound) %>% summarize(count = n(), coloc_variants = paste(locus_tmp, collapse = '; ')) %>% left_join(., compound_info, by = c('compound' = 'id')) %>% arrange(type, desc(count)) %>% ungroup

# How many loci colocalise with a metabolite
compound_count_loci %>% filter(type == 'plasma') %>% pull(count) %>% summary
compound_count_loci %>% filter(type == 'urine') %>% pull(count) %>% summary

# How many compounds colocalise at a locus
loci_count_met = all_sig %>% mutate(locus_tmp = gsub('_[A-Z]', '', locus)) %>% distinct(compound, locus_tmp, type) %>% group_by(type, locus_tmp) %>% summarize(count = n(), compound_list = paste(compound, collapse = '; ')) %>% arrange(type, desc(count)) %>% ungroup

loci_count_met %>% filter(type == 'plasma') %>% pull(count) %>% summary
loci_count_met %>% filter(type == 'urine') %>% pull(count) %>% summary

top_plasma = all_sig %>% filter(type == 'plasma', locus == 'chr2_26.91_28.71MB')
top_urine = all_sig %>% filter(type == 'urine', grepl('chr6_25.07_32.85MB', locus))

# What is the overlap between plasma and urine (i.e. which compound showed up in both urine and plasma)
compound_count_type = all_sig %>% distinct(type, phenostring) %>% group_by(phenostring) %>% summarize(count = n(), type_list = paste(sort(type), collapse = '/'))

table(compound_count_type$type_list)

both = compound_count_type %>% filter(type_list == 'plasma/urine')
plasma_only = compound_count_type %>% filter(type_list == 'plasma')
urine_only = compound_count_type %>% filter(type_list == 'urine')

# How many gout loci did urate colocalise with
compound_count_loci %>% filter(phenostring == 'urate') %>% as.data.frame

table(compound_count_loci$type, compound_count_loci$count)

# Save a bunch of files
write.table(full_res, 'results/merged_res/coloc_merged.full.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(all_sig, 'results/merged_res/coloc_merged.full.pp0.8.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(compound_count_var, 'results/merged_res/snp_per_met.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(compound_count_loci, 'results/merged_res/loci_per_met.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(loci_count_met, 'results/merged_res/met_per_loci.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(both, 'results/merged_res/plasma_urine_overlap.pp0.8.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(plasma_only, 'results/merged_res/plasma_only.pp0.8.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(urine_only, 'results/merged_res/urine_only.pp0.8.txt', sep = '\t', col.names = T, row.names = F, quote = F)

