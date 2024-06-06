# Quick script to narrow down which metabolites to use for the MR analysis
#
# - Take the 144 plasma metabolites
# - Group by metabolites and count up loci
# - Use/select metabolites that have X colocalised loci or more

library(dplyr)
library(tidyr)
library(purrr)

setwd('..')

# Load in list of 144 common plasma metabolites
met = readLines('results/met_overlap/all_plasma_overlap.txt')

# Load all the plasma metabolites:
metsim = read.table('results/coloc_res/metsim/merged_res/coloc_merged.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
metsim_hmdb = read.table('results/hmdb_match/metsim/all_metsim.hmdb_merge.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
metsim = left_join(metsim, metsim_hmdb, by = c('compound' = 'phenocode', 'category', 'phenostring'), suffix = c('', '.hmdb'), relationship = 'many-to-many')

schl = read.table('results/hmdb_match/schlosser/coloc_results.hmdb_merge.txt', sep = '\t', header = T, stringsAsFactors = F)
schl = schl %>% filter(category == 'plasma')
chen = read.table('results/hmdb_match/chen/coloc_results.hmdb_merge.txt', sep = '\t', header = T, stringsAsFactors = F)

dat_list = list(metsim, schl, chen)

# Load and merge the master phenostring information so I can pull out the 144
# metabolites
master = read.table('results/hmdb_match/master_phenostring.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

dat_list = map(dat_list, ~ left_join(.x, master, by = 'metabolite.hmdb', suffix = c('.original', '')))

# Pull out the 144 plasma metabolites
dat_sub = map(dat_list, ~ .x %>% filter(metabolite %in% met))

# For Megan:
# megan = map_dfr(dat_list, ~ .x %>% select(locus, SNP, metabolite)) %>% distinct
# write.table(megan, 'results/met_overlap/met_SNP_pairs.for_megan.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Group by metabolite, then count up the number of unique SNP/loci that
# colocalised with each metabolite
dat_count = map(dat_sub, ~ .x %>% distinct(SNP, metabolite.original) %>% group_by(metabolite.original) %>% summarise(coloc_loci = n()) %>% arrange(desc(coloc_loci)))

res = map(dat_count, ~ .x %>% filter(coloc_loci >= 4) %>% pull(metabolite.original))

# Pull out the HMDB metabolite name so I can pull out the study-specific
# metabolite ID (GCST*, C*, etc.) for all the metabolites
res_hmdb = map2(dat_sub, res, ~ .x %>% filter(metabolite.original %in% .y) %>% pull(metabolite.hmdb))
res_hmdb = unique(unlist(res_hmdb))

# Write out the list of metabolites that had at least 4 colocalising loci with
# gout.
# Probably best to pull out the compound ID, so I can use the list to directly
# refer to the metQTL filename

out_names = 'results/mr_results/STUDY.mr_list.txt'
out_names = map_chr(c('metsim', 'schl', 'chen'), ~ gsub('STUDY', .x, out_names))

# NOTE: 2-hydroxybutyrate and 2-hydroxyisobutyrate are the same compound
tmp_res = map(list(metsim, schl, chen), ~ .x %>% filter(metabolite.hmdb %in% res_hmdb) %>% pull(compound) %>% unique)

map2(tmp_res, out_names, ~ writeLines(.x, .y))

# Now save the study-specific metabolites:
tmp_res = map2(list(metsim, schl, chen), res,  ~ .x %>% filter(metabolite %in% .y) %>% pull(compound) %>% unique)

out_names = map(out_names, ~ gsub('.txt', '.specific.txt', .x))

# map2(tmp_res, out_names, ~ writeLines(.x, .y))

# tmp = metsim %>% select(metabolite.hmdb, compound) %>% distinct
# tmp = schl %>% select(metabolite.hmdb, compound) %>% distinct %>% left_join(tmp, ., by = 'metabolite.hmdb', suffix = c('', '.schl'))
# tmp = chen %>% select(metabolite.hmdb, compound) %>% distinct %>% left_join(tmp, ., by = 'metabolite.hmdb', suffix = c('.metsim', '.chen'))
# tmp = tmp %>% filter(metabolite.hmdb %in% res_hmdb)

# tmp %>% filter(compound.chen %in% c("GCST90199621", "GCST90199971", "GCST90200307", "GCST90199666", "GCST90200065"))
