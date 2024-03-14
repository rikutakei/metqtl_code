# NOTE: the HMDB version is version 5, released 2021-11-17
library(dplyr)
library(purrr)
library(tidyr)
library(fedmatch)

################################################################################
# Libraries for MetaboAnalystR (ended up not using it, but here for the record)

# library(BiocManager)
# library(devtools)
# pkgs = c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "crmn", "httr", "qs")
# BiocManager::install(pkgs)
# devtools::install_github('xia-lab/MetaboAnalystR', build = T, build_vignettes = T, build_manual = T)
# library(MetaboAnalystR)

################################################################################

# setwd('..')

# Load colocalisation results from gout vs metQTL
dat = read.table('results/coloc_res/metsim/merged_res/coloc_merged.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

# First remove the "Uncharacterized" metabolites
dat = dat %>% filter(!(category %in% c('Uncharacterized', 'Partially_Characterized')))

# List to use on the MetaboAnalyst website to see if matches are better there
# writeLines(unique(dat$phenostring[dat$category != 'Lipid']), 'results/coloc_res/metsim/merged_res/metaboanalyst_list.met.txt')
# writeLines(unique(dat$phenostring[dat$category == 'Lipid']), 'results/coloc_res/metsim/merged_res/metaboanalyst_list.lipid.txt')

# Need to clean the metabolite name a little
# - remove "*" character
# - remove the secondary name that is in bracket (e.g. "phenyllactate (PLA)")
# - separate out the alternative metabolite after the slash (e.g. "oleate/vaccenate")
dat$metabolite = gsub('\\*.*', '', dat$phenostring)
dat$metabolite = gsub(' \\[.\\]', '', dat$metabolite)
dat$metabolite = gsub(' \\([^\\)]*\\) *$', '', dat$metabolite)
dat$metabolite = gsub(' \\(.*\\)$', '', dat$metabolite)
dat = dat %>% separate_rows(metabolite, sep = '/') %>% as.data.frame

# Aconitate has a "[cis or trans]" in the phenostring - separate this into cis-
# and trans-aconitate
ind = which(grepl('cis or trans', dat$metabolite))

if (length(ind) > 0) {
  tmp = dat[c(ind, ind), ]
  tmp$metabolite[1] = 'cis-aconitate'
  tmp$metabolite[2] = 'trans-aconitate'
  dat = dat[-c(ind), ]
  dat = rbind(dat, tmp)
}

# Lower case all the metabolites
dat$metabolite = tolower(dat$metabolite)

# Load in HMDB data for matching the metQTL metabolites with HMDB
met_dat = read.table('data/hmdb/hmdb_synonyms.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_dat$metabolite = tolower(met_dat$metabolite)
met_dat$synonym = tolower(met_dat$synonym)

# Add on "status" info:
status = read.table('data/hmdb/hmdb_status.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
status$metabolite = tolower(status$metabolite)
status$status_value = map_dbl(status$status, ~ switch(.x, 'quantified' = 1, 'detected' = 2, 'expected' = 3, 'predicted' = 4, 100))

met_dat = left_join(met_dat, status) %>% distinct

# Only need to match the unique metabolites in the coloc result:
sub_dat = dat %>% distinct(metabolite)

# Merge the two data sets for exact matching, and then use fuzzy matching for
# those that didn't have exact matches

# NOTE: Since some amino acids have D- and L- form (or just duplicate entry),
# keep both in the merged data set since different isoforms have different
# pathways/proteins associated with it. With that said, there are clearly some
# wrong metabolites in the synonym list in HMDB.
exact = left_join(sub_dat, met_dat, by = c('metabolite' = 'synonym'), suffix = c('', '.hmdb'), relationship = 'many-to-many', keep = T) %>% group_by(metabolite) %>% arrange(status_value) %>% ungroup
exact = exact %>% filter(!is.na(metabolite.hmdb)) %>% as.data.frame

# Due to the "simplified" metabolite name in my data set and multiple
# metabolite matches from the synonyms in HMDB, I need to weed out the
# irrelevant metabolite entries from HMDB (especially the lipids, like
# diacylglecerol)
dup_pheno = exact$metabolite[which(duplicated(exact$metabolite))]
dup = exact %>% filter(metabolite %in% dup_pheno) %>% arrange(metabolite)
non_dup = exact %>% filter(!(metabolite %in% dup_pheno))

# Good starting point is to filter out by the "status" of the metabolite in
# HMDB - the "quantified" metabolites tend to have more complete information,
# then "detected", "expected", and lastly "predicted".
#
# Find and keep the metabolite(s) with best status for each query metabolite
# and keep ties
target = dup %>% group_by(metabolite) %>% arrange(status_value) %>% slice(1) %>% ungroup %>% distinct(metabolite, status_value)

dup_clean = inner_join(dup, target)

# There can be metabolites with candidate matches that have the same status.
# These tend to be L-/D- form of the same compound (so only two candidates are
# present), or a large number of candidates due to being lipids.
#
# Filter out those that have more than 2 candidates:
dup_met = dup_clean %>% group_by(metabolite) %>% summarise(count = n()) %>% filter(count > 2) %>% pull(metabolite)
good = dup_clean %>% filter(!(metabolite %in% dup_met))

dup_clean = dup_clean %>% filter(metabolite %in% dup_met)

# Now, check for exact match between the metabolite names from my data vs HMDB
# and add it to the good list
name_match = dup_clean %>% filter(metabolite == metabolite.hmdb)

if (nrow(good) > 0) {
  good$match.type = 'exact_match'
}

if (nrow(name_match) > 0) {
  name_match$match.type = 'exact_match'
}

res = rbind(good, name_match)

dup_clean = dup_clean %>% filter(!(metabolite %in% res$metabolite))

# Now I'm left with the metabolites that don't quite match with the HMDB
# metabolite and also have the same synonyms across multiple metabolites.
#
# Try fuzzy matching to reduce the table into candidates that have similar
# names (e.g. D-/L-amino acids)
dat1 = dup_clean %>% select(!metabolite)
dat2 = dup_clean %>% select(!metabolite.hmdb) %>% distinct
dat1$id1 = rownames(dat1)
dat2$id2 = rownames(dat2)

fuzzy = merge_plus(data1 = dat1, data2 = dat2, by.x = 'metabolite.hmdb', by.y = 'metabolite', match_type = 'fuzzy', unique_key_1 = 'id1', unique_key_2 = 'id2', fuzzy_settings = build_fuzzy_settings(p = 0.25, maxDist = 0.25))

# Only interested in the matches that have the same synonyms:
fuzzy = fuzzy$matches %>% filter(synonym_1 == synonym_2) %>% arrange(metabolite)

# Pull out metabolites that got narrowed down into a single match and add it to
# the result:
dup_met = fuzzy$metabolite[which(duplicated(fuzzy$metabolite))]
fuz_match = fuzzy %>% filter(!(metabolite %in% dup_met)) %>% select(metabolite, metabolite.hmdb, synonym_1, status_1, status_value_1)
colnames(fuz_match) = gsub('_1', '', colnames(fuz_match))

fuz_match$match.type = 'fuzzy_match'

res = rbind(res, fuz_match)

# For the other metabolites that had fuzzy matches, the duplicate entries may
# still be valid (e.g. D-/L-amino acids) - save and check manually
fuz_check = fuzzy %>% filter(metabolite %in% dup_met) %>% select(metabolite, metabolite.hmdb, synonym_1)

write.table(fuz_check, 'results/hmdb_match/fuzzy_match_check.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Read in the metabolites to keep/add to `res`
fuz_keep = read.table('results/hmdb_match/fuzzy_match_keep.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
fuz_keep = inner_join(fuzzy, fuz_keep) %>% select(metabolite, metabolite.hmdb, synonym_1, status_1, status_value_1)
colnames(fuz_keep) = gsub('_1', '', colnames(fuz_keep))

fuz_keep$match.type = 'fuzzy_match'

res = rbind(res, fuz_keep)

# For the rest of the metabolites that didn't match well through fuzzy matching
# or didn't pass the manual checking, write it out and check manually (these
# are mostly lipids). Also, ignore the "status" filter that was used earlier,
# since some of the better matches were included in the "lower" status and got
# filtered out
dup_clean = dup_clean %>% filter(!(metabolite %in% res$metabolite))
dup_met = unique(dup_clean$metabolite)
check = exact %>% filter(metabolite %in% dup_met)
check = dat %>% distinct(phenostring, metabolite) %>% left_join(., check, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb)) %>% arrange(metabolite)

write.table(check, 'results/hmdb_match/lipid_multi_check.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Read in the metabolites to keep:
keep = read.table('results/hmdb_match/lipid_multi_keep.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
keep$match.type = 'manual'

# Add phenostring back to `res` and rbind the `keep`
res = dat %>% distinct(phenostring, metabolite) %>% left_join(res, ., relationship = 'many-to-many') %>% select(phenostring, metabolite:match.type)
res = rbind(res, keep)

# Now merge the cleaned up duplicate metabolite table with the rest of
# non-duplicated results:
non_dup = dat %>% distinct(phenostring, metabolite) %>% left_join(non_dup, ., relationship = 'many-to-many') %>% select(phenostring, metabolite:status_value)
non_dup$match.type = 'exact_match'

res = rbind(res, non_dup)

################################################################################
# Now focus on the metabolites that didn't match exactly with the synonyms (or
# did match, but none were correct - e.g. some diacylglycerols)

no_match = sub_dat %>% filter(!(metabolite %in% res$metabolite))

# Need a unique ID in both data sets for fuzzy matching:
met_dat$id = rownames(met_dat)
no_match$uniq_id = rownames(no_match)

fuzzy_res = merge_plus(met_dat, no_match, by.x = 'synonym', by.y = 'metabolite', match_type = 'fuzzy', unique_key_1 = 'id', unique_key_2 = 'uniq_id')

# Pull out the `no_match` metabolites:
cand = fuzzy_res$matches %>% filter(metabolite.y %in% no_match$metabolite) %>% arrange(as.numeric(uniq_id))
colnames(cand)[c(3, 7)] = c('metabolite.hmdb', 'metabolite')

# Write it out and check manually:
write.table(cand, 'results/hmdb_match/non_match_candidates.check.txt', sep = '\t', col.names = T, row.names = F, quote = F)

cand_keep = read.table('results/hmdb_match/non_match_candidates.keep.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

cand_keep = dat %>% distinct(phenostring, metabolite) %>% left_join(cand_keep, ., relationship = 'many-to-many') %>% select(phenostring, metabolite, metabolite.hmdb, synonym, status, status_value)
cand_keep$match.type = 'manual'

res = rbind(res, cand_keep)

# Write out the metabolites that didn't match well/good enough with the HMDB:
no_match = sub_dat %>% filter(!(metabolite %in% res$metabolite))
tmp = unique(no_match$metabolite)
tmp = tmp[!grepl('x-', tmp)]
writeLines(tmp, 'results/hmdb_match/no_match_met.txt')

dat = left_join(dat, res, relationship = 'many-to-many')

clean = dat %>% filter(!is.na(match.type))

write.table(dat, 'results/hmdb_match/coloc_results.hmdb_merge.txt', sep = '\t', col.names = T, row.names = F)
write.table(clean, 'results/hmdb_match/coloc_results.hmdb_merge.match_only.txt', sep = '\t', col.names = T, row.names = F)

