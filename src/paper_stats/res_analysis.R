library(dplyr)
library(tidyr)

setwd('..')

compound_info = read.delim('data/metqtl_pheno/full_pheno.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
compound_info = compound_info %>% select(phenostring:summaryStatistics.urine) %>% pivot_longer(cols = c(phenocode:summaryStatistics.urine), names_to = 'dummy', values_to = 'id') %>% filter(id != '0') %>% select(phenostring, id) %>% distinct

full_res = read.table('results/merged_res/coloc_merged.full.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
all_sig = read.table('results/merged_res/coloc_merged.full.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
loci_count_met = read.table('results/merged_res/met_per_loci.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
both = read.table('results/merged_res/plasma_urine_overlap.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
urine_only = read.table('results/merged_res/urine_only.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
compound_count_var = read.table('results/merged_res/snp_per_met.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
compound_count_loci = read.table('results/merged_res/loci_per_met.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
plasma_only = read.table('results/merged_res/plasma_only.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

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

# There are 10 compounds that showed coloc with at least 5 gout loci in either
# plasma or urine:
cand = compound_count_loci %>% filter(count >= 5)
cand %>% pull(phenostring) %>% sort

# How many compounds showed coloc with at least 5 gout loci if plasma and urine
# metabolites are considered together:
#
# NOTE: Since the `id` is different between plasma and urine, you need to group
# by the `phenostring` column
top_comp = all_sig %>% mutate(locus_tmp = gsub('_[A-Z]', '', locus)) %>% distinct(phenostring, locus_tmp) %>% group_by(phenostring) %>% summarize(count = n(), coloc_variants = paste(locus_tmp, collapse = '; ')) %>% arrange(desc(count)) %>% ungroup

table(top_comp$count)

# Not separating the plasma/urine metabolites makes more sense for what I am
# doing ("what genetic loci affects the metabolite") - if it was digging deeper
# about which locus affecting plasma and/or urine side of things, then it would
# probably matter
#
# I arbitrarily chose 5 loci as a cutoff before, but it would make more sense
# if I used urate as the baseline/reference point, since we already know urate
# is causal of gout. Urate colocalises with 4 gout loci, so use that as the
# cutoff:
cand = top_comp %>% filter(count >= 4)

# Pull out the compound ID (for both plasma and urine):
compound_info = read.delim('data/metqtl_pheno/full_pheno.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
cand_id = compound_info %>% filter(phenostring %in% cand$phenostring)
cand_id$plasma_id = paste(cand_id$phenocode, cand_id$summaryStatistics.plasma, sep = '_')

urine_files = cand_id %>% filter(summaryStatistics.urine != '0') %>% pull(summaryStatistics.urine)
urine_files = paste('data/schlosser_metqtl/urine/', urine_files, sep = '')
urine_files = paste(urine_files, '_buildGRCh37.tsv.gz', sep = '')

plasma_metafiles = cand_id %>% filter(!grepl('^0_|_0$', plasma_id)) %>% pull(plasma_id)
plasma_metafiles = paste('data/plasma_meta/', plasma_metafiles, sep = '')
plasma_metafiles = paste(plasma_metafiles, '.meta1.tbl.gz', sep = '')

plasma_files = cand_id %>% filter(grepl('^0_|_0$', plasma_id)) %>% pull(plasma_id)
plasma_files = gsub('^0_|_0$', '', plasma_files)
ind1 = which(grepl('^G', plasma_files))
ind2 = which(grepl('^C', plasma_files))
plasma_files[ind1] = paste('data/schlosser_metqtl/plasma/', plasma_files[ind1], sep = '')
plasma_files[ind1] = paste(plasma_files[ind1], '_buildGRCh37.tsv.gz', sep = '')
plasma_files[ind2] = paste('data/metsim_data/', plasma_files[ind2], sep = '')
plasma_files[ind2] = paste(plasma_files[ind2], '.clean.txt.gz', sep = '')

cand_files = c(urine_files, plasma_metafiles, plasma_files)

writeLines(cand_files, 'data/mr/cand_compound_for_mr.txt')

