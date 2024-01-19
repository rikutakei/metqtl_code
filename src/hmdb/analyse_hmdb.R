# Analyze colocalised metabolites that have matching HMDB metabolites
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)

setwd('../..')

# Load in HMDB-matched coloc results:
dat = read.table('results/hmdb_match/coloc_results.hmdb_merge.match_only.txt', sep = '\t', header = T, stringsAsFactors = F)

# Take a look at the classes of metabolites:
met_class = read.table('data/hmdb/hmdb_chemtax.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_class = met_class %>% mutate(metabolite = tolower(metabolite))

# Add to dat and clean it up for plotting
dat_class = left_join(dat, met_class, by = c('metabolite.hmdb' = 'metabolite')) %>% select(metabolite, type, kingdom:molecular_framework)

# Since a lot of metabolites are in common between plasma and urine, make
# a third type for both
dat_class$type = ifelse(dat_class$type == 'urine', 'Urine-only', 'Plasma-only')
dat_class = dat_class %>% group_by(metabolite) %>% mutate(type = paste(unique(type), collapse = '/')) %>% ungroup
dat_class$type[which(grepl('/', dat_class$type))] = 'Both'
dat_class = dat_class %>% distinct

dat_class = dat_class %>% pivot_longer(cols = kingdom:molecular_framework, names_to = 'category', values_to = 'value') %>% filter(!is.na(value), value != '') %>% distinct

# Kingdom is useless, as they are all organic compounds:
dat_class %>% filter(category == 'kingdom') %>% ggplot(aes(x = value)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))

dat_class %>% filter(category == 'super_class') %>% ggplot(aes(x = value, fill = type)) +
  geom_bar() +
  facet_grid(type ~ .) +
  theme_bw() +
  labs(x = 'Super class', y = 'Count') +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(hjust = 0.95, vjust = 0.5, angle = 90)
  )

ggsave('results/plots/hmdb_superclass.png', width = 6, height = 12)

dat_class %>% filter(category == 'class') %>% ggplot(aes(x = value, fill = type)) +
  geom_bar() +
  facet_grid(type ~ .) +
  theme_bw() +
  labs(x = 'Class', y = 'Count') +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(size = 5, hjust = 0.95, vjust = 0.5, angle = 90)
  )

ggsave('results/plots/hmdb_class.png', width = 6, height = 12)

dat_class %>% filter(category == 'sub_class') %>% ggplot(aes(x = value, fill = type)) +
  geom_bar() +
  facet_grid(type ~ .) +
  theme_bw() +
  labs(x = 'Sub class', y = 'Count') +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(size = 5, hjust = 0.95, vjust = 0.5, angle = 90)
  )

ggsave('results/plots/hmdb_subclass.png', width = 6, height = 12)

dat_class %>% filter(category == 'molecular_framework') %>% ggplot(aes(x = value, fill = type)) +
  geom_bar() +
  facet_grid(type ~ .) +
  theme_bw() +
  labs(x = 'Direct parent', y = 'Count') +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(hjust = 0.95, vjust = 0.5, angle = 90)
  )

ggsave('results/plots/hmdb_parent.png', width = 6, height = 12)

# Take a look at which proteins/genes are involved with each metabolite:
met_prot = read.table('data/hmdb/hmdb_metprot.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_prot = met_prot %>% mutate(metabolite = tolower(metabolite))

dat_prot = left_join(dat, met_prot, by = c('metabolite.hmdb' = 'metabolite'), relationship = 'many-to-many') %>% filter(!is.na(uniprot_id)) %>% distinct(locus, SNP, type, phenostring, metabolite.hmdb, name, uniprot_id, gene_name, protein_type)

# Load in lead SNP info
loci = read.table('data/gwas/indep_snps_from_paper.clean.txt', sep = '\t', header = T, stringsAsFactors = F)

dat_prot = loci %>% select(!cohort) %>% distinct %>% left_join(dat_prot, .) %>% mutate(start = pos - 500000, end = pos + 500000)

gene_info = read.table('data/gene_info/UCSC_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactors = F)
gene_info = gene_info %>% filter(!(Chrom %in% c('chrM', 'chrY')), !grepl('_', Chrom)) %>% select(Chrom:End, Gene) %>% distinct
gene_info = gene_info %>% group_by(Gene) %>% slice(1) %>% ungroup

dat_prot = left_join(dat_prot, gene_info, by = c('gene_name' = 'Gene'))

# Now filter out genes that are not within (or overlap with) the coloc region.
# First, those not on the same chromosome
dat_prot = dat_prot %>% mutate(chr = paste('chr', chr, sep = ''))
dat_prot$chr = gsub('23', 'X', dat_prot$chr)
dat_prot = dat_prot %>% filter(chr == Chrom)

# Check positioning of the gene.
# If the SNP position, locus start, or locus end is within the gene region,
# then the gene (or at least part of the gene) is within the locus
dat_prot = dat_prot %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End))

# Save for now:
write.table(dat_prot, 'results/hmdb_match/coloc_results.met_genes.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Analyses that can be done:
# - Pull out genes
#   - check if genes within locus
#   - check if genes encode proteins in the same pathway
# - Check if metabolites from a single locus is part of one (or more) biological pathway(s)

# Majority of the colocalised metabolites are NOT directly
# controlled/transported by the protein product of the gene at the locus.
# It is likely that the metabolite is either a stable intermediate or an
# end-product of a pathway, where the gene(s) at the locus affects the
# metabolites before and/or after the colocalised metabolite. Also, HMDB is
# a little "patchy" in terms of transporters (e.g. doesn't have ABCG2 or
# SLC22A9 listed for uric acid)
#
# See if any of the proteins/genes involved in any of the metabolites that is
# part of any of the pathways with the colocalised metabolites is within the
# locus.

# Load in metabolite-pathway data:
met_path = read.table('data/hmdb/hmdb_biopath.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_path = met_path %>% mutate(metabolite = tolower(metabolite))

# Load in metabolite-class data:
met_class = read.table('data/hmdb/hmdb_chemtax.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_class = met_class %>% mutate(metabolite = tolower(metabolite))

# Pull out all the metabolites from the colocalisation results:
met_list = dat %>% pull(metabolite.hmdb) %>% unique

# Pull out all the pathways that have one of the metabolites in the list
path_list = met_path %>% filter(metabolite %in% met_list) %>% pull(name)

# Now pull out all of the metabolites that are involved in any of the pathways
# listed above
full_list = met_path %>% filter(name %in% path_list)

# Format the table so that each metabolite has a pathway entry, and then a list
# of metabolites involved in the listed pathway:
full_list = full_list %>% select(metabolite, name) %>% left_join(full_list, ., by = 'name', suffix = c('', '.pathway'), relationship = 'many-to-many') %>% distinct

# Exclude "metabolites" that are inroganic:
inorg = met_class %>% filter(kingdom == 'Inorganic compounds') %>% select(metabolite, kingdom)
full_list = left_join(full_list, inorg, by = 'metabolite')
full_list = left_join(full_list, inorg, by = c('metabolite.pathway' = 'metabolite'), suffix = c('.met', '.pth'))
full_list = full_list %>% filter(is.na(kingdom.met) & is.na(kingdom.pth)) %>% select(!contains('kingdom'))

# Remove metabolites that are present in the pathways of >= X% of metabolites.
#
# For example, ATP ("adenosine triphosphate") is a part of at least one pathway
# in 91.8% of the metabolites - the odds of linking ATP with a given metabolite
# is high, so remove them to reduce noise in the end result
perc_met = full_list %>% group_by(metabolite) %>% summarize(total = length(unique(metabolite.pathway))) %>% ungroup %>% mutate(perc_met = total / nrow(.))
full_list = left_join(full_list, perc_met, by = 'metabolite')

# There is a very sharp drop-off in the percentage after 15%, where the next
# percentage after that is 55%:
perc_met %>% arrange(desc(perc_met)) %>% filter(perc_met >= 0.03) %>% distinct(metabolite, perc_met) %>% as.data.frame
full_list %>% ggplot(aes(perc_met)) +
  geom_density()

# Remove metabolites that are part of >10% of other metabolites' pathway:
hi_perc = full_list %>% filter(perc_met >= 0.1) %>% pull(metabolite) %>% unique
full_list = full_list %>% filter(!(metabolite %in% hi_perc) & !(metabolite.pathway %in% hi_perc))

# Reduce the list down to those present in the data
full_list = full_list %>% filter(metabolite %in% dat$metabolite.hmdb) %>% distinct

# Now link up the enzymes/proteins/genes involved with each metabolite:
full_list = left_join(full_list, met_prot, by = c('metabolite.pathway' = 'metabolite'), suffix = c('.pathway', '.protein'), relationship = 'many-to-many')

# Make a smaller table for merging back the pathways based on metabolites and
# proteins:
path_prot = full_list %>% select(metabolite, name.pathway, metabolite.pathway, uniprot_id, gene_name) %>% distinct

# Need to make the `full_list` smaller in order to merge with `dat` - this is
# done by ignoring the pathway column in which the metabolite is associated
# with (hence `path_prot` was generated beforehand to merge it back in):
full_list = full_list %>% filter(!is.na(gene_name)) %>% select(metabolite, metabolite.pathway:protein_type) %>% distinct

# Repeat the process for the metabolite proteins, but now with `full_list`
dat_prot = left_join(dat, full_list, by = c('metabolite.hmdb' = 'metabolite'), relationship = 'many-to-many') %>% filter(!is.na(gene_name))

dat_prot = dat_prot %>% distinct(locus, SNP, type, phenostring, metabolite.hmdb, metabolite.pathway, name.protein, uniprot_id, gene_name, protein_type)

dat_prot = loci %>% select(!cohort) %>% distinct %>% left_join(dat_prot, .) %>% mutate(start = pos - 500000, end = pos + 500000)

dat_prot = left_join(dat_prot, gene_info, by = c('gene_name' = 'Gene'))

# Now filter out genes that are not within (or overlap with) the coloc region.
# First, those not on the same chromosome
dat_prot = dat_prot %>% mutate(chr = paste('chr', chr, sep = ''))
dat_prot$chr = gsub('23', 'X', dat_prot$chr)
dat_prot = dat_prot %>% filter(chr == Chrom)

# Check positioning of the gene.
# If the SNP position, locus start, or locus end is within the gene region,
# then the gene (or at least part of the gene) is within the locus
dat_prot = dat_prot %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End)) %>% as.data.frame

# How many genes at how many loci?
length(unique(dat_prot$gene_name))
length(unique(dat_prot$locus))

write.table(dat_prot, 'results/hmdb_match/coloc_results.met_pathway_genes.txt', sep = '\t', col.names = T, row.names = F, quote = F)

dat_prot$type = ifelse(dat_prot$type == 'urine', 'Urine-only', 'Plasma-only')
dat_prot = dat_prot %>% group_by(metabolite.hmdb) %>% mutate(type = paste(unique(type), collapse = '/')) %>% ungroup
dat_prot$type[which(grepl('/', dat_prot$type))] = 'Both'

# How many different genes/proteins (that are also within the colocalised
# region) act on the pathway that involves the colocalised metabolite
dat_prot %>% filter(!grepl('dg', metabolite.pathway)) %>%
  distinct(phenostring, type, gene_name) %>% ggplot(aes(x = str_to_title(phenostring), fill = type)) +
  geom_bar() +
  facet_grid(type ~ .) +
  theme_bw() +
  labs(x = 'Metabolite', y = 'Number of proteins in a pathway involving the metabolite that are also present at the colocalised region') +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  )

ggsave('results/plots/hmdb_met.gene_count.png', width = 16, height = 12)

# Need to add some stats to the genes that show up at the colocalised locus.
#
# Draw up a 2-by-2 contingency table, where the categories are:
# - Whether a gene is present at the colocalised locus; and
# - Whether a gene is part of the colocalised metabolite's pathway
#
# Total number will bet the genes present in the HMDB

# First, need to add gene start/end positions to determine how many genes are
# present at each locus, and to also determine how many total (usable) genes
# are in HMDB:
hmdb_gene = met_prot %>% distinct(gene_name) %>% left_join(., gene_info, by = c('gene_name' = 'Gene')) %>% filter(!is.na(Chrom))

# How many genes in total? -> 5594 genes
total_genes = length(unique(hmdb_gene$gene_name))

# Figure out what and how many genes overlap with each locus:
locus_gene = dat_prot %>% select(locus, chr:end) %>% distinct
locus_gene = left_join(locus_gene, hmdb_gene, by = c('chr' = 'Chrom'), relationship = 'many-to-many')
locus_gene = locus_gene %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End)) %>% as.data.frame
locus_gene_count = locus_gene %>% distinct(locus, gene_name) %>% group_by(locus) %>% summarize(total = n()) %>% ungroup %>% distinct

# Need a table of colocalised metabolite with all the genes in the associated
# pathways. Use the `full_list` data which has been filtered down to the
# metabolites that colocalised:
coloc_met_path = full_list %>% distinct(metabolite, gene_name)

# Function to generate a 2-by-2 contingency table:
make_cont = function(locus_name, metabolite_name, gene, total = 5594) {
  # head(dat_prot) %>% as.data.frame
  # Pull out genes in locus
  tmp1 = locus_gene %>% filter(locus == locus_name) %>% pull(gene_name)
  # Pull out genes from metabolite's pathways:
  tmp2 = coloc_met_path %>% filter(metabolite == metabolite_name)
  # How many genes in the locus are in the metabolite's pathway?
  tmp3 = which(tmp1 %in% tmp2$gene_name)
  # Make contingency table:
  a = length(tmp3)
  b = length(tmp1) - a
  c = nrow(tmp2) - a
  d = total - sum(c(a, b, c))
  tab = matrix(c(a, b, c, d), ncol = 2)
  return(tab)
}

cont_tab = pmap(list(dat_prot$locus, dat_prot$metabolite.hmdb, dat_prot$gene_name), ~ make_cont(..1, ..2, ..3, total_genes))
# test = pmap(list(head(dat_prot$locus), head(dat_prot$metabolite.hmdb), head(dat_prot$gene_name)), ~ make_cont(..1, ..2, ..3, total_genes))
# test = pmap(list(head(dat_prot$locus), head(dat_prot$metabolite.hmdb), head(dat_prot$gene_name)), ~ head)

head(cont_tab)

# Run Fisher's exact test on the tables:
fish_test = map_dbl(cont_tab, ~ fisher.test(.x)$p.value)
dat_prot$p.fisher = fish_test
dat_prot$p.fdr = p.adjust(dat_prot$p.fisher, method = 'fdr')

# Load in the MR results to filter out certain metabolites
mr = read.table('results/mr_results/mr_res.p_summary.txt', sep = '\t', header = T, stringsAsFactors = F)

# Make a list of all significant MR metabolites (both metQTL and gout) and
# a list with just those that were significant with no evidence of pleiotropy
# via MR-Egger
mr_sig = mr %>% filter(p.IVW <= 0.05 | p.WM <= 0.05)
mr_sig_nop = mr_sig %>% filter(p.egger_int > 0.05)

# Pull out significantly enriched pathways
dat_prot_mr = dat_prot %>% filter(phenostring %in% mr_sig$phenostring)
dat_prot_mr_nop = dat_prot %>% filter(phenostring %in% mr_sig_nop$phenostring)

dat_prot_mr %>% as.data.frame %>% head

# TODO: from here
unique_test = 
threshold = 0.05 / nrow(dat_prot_mr)
sig_path = dat_prot_mr %>% filter(p.fisher <= threshold) %>% pull(gene_name) %>% unique
sig_path_fdr = dat_prot_mr %>% filter(p.fdr <= 0.05) %>% pull(gene_name) %>% unique
sig_path_fdr

threshold = 0.05 / nrow(dat_prot_mr_nop)
sig_path_nop = dat_prot_mr_nop %>% filter(p.fisher <= threshold) %>% pull(gene_name) %>% unique
sig_path_fdr_nop = dat_prot_mr_nop %>% filter(p.fdr <= 0.05) %>% pull(gene_name) %>% unique

threshold = 0.05 / nrow(dat_prot)
dat_prot %>% filter(p.fisher <= 0.05) %>% pull(gene_name) %>% unique
dat_prot %>% filter(p.fisher <= 0.05) %>% pull(locus) %>% unique
dat_prot %>% filter(p.fisher <= 0.05) %>% pull(metabolite.hmdb) %>% unique
dat_prot %>% filter(p.fdr <= 0.05) %>% pull(gene_name) %>% unique
dat_prot %>% filter(p.fdr <= 0.05) %>% pull(locus) %>% unique
dat_prot %>% filter(p.fdr <= 0.05) %>% pull(metabolite.hmdb) %>% unique
dat_prot %>% filter(p.fisher <= threshold) %>% pull(gene_name) %>% unique
dat_prot %>% filter(p.fisher <= threshold) %>% pull(locus) %>% unique
dat_prot %>% filter(p.fisher <= threshold) %>% pull(metabolite.hmdb) %>% unique

################################################################################

# What "metabolite" shows up the most when you consider the whole pathway
dat_prot %>% filter(!grepl('dg', metabolite.pathway)) %>%
  distinct(metabolite.pathway, type, gene_name) %>% ggplot(aes(x = str_to_title(metabolite.pathway), fill = type)) +
  geom_bar() +
  facet_grid(type ~ .) +
  theme_bw() +
  labs(x = 'Metabolite', y = 'Count') +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  )

ggsave('results/plots/hmdb_met.path_met.png', width = 16, height = 12)

# TODO: 
# Merge back in the pathway(s) in which the metabolite came from, but first

tmp = left_join(dat_prot, path_prot, by = c('metabolite.hmdb' = 'metabolite', 'metabolite.pathway', 'gene_name', 'uniprot_id'), relationship = 'many-to-many') %>% distinct
tmp = tmp %>% select(locus:name.protein, gene_name, name.pathway)

dim(dat_prot)
dim(tmp)
head(tmp) %>% as.data.frame


write.table(tmp, 'results/hmdb_match/coloc_results.met_pathway_genes.tmp.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# TODO: The problem with the above approach is that you get way too many
# pathways, even if you narrow it down to matching metabolite and gene/protein
# and it is also hard to figure out what pathway the metabolite is part of.
# Instead, try and pull out the metabolite(s) the gene/protein is involved
# with, pull out the pathways for each metabolite, and then see what pathway(s)
# are the most common among the metabolites that the gene is involved with.

# Function to pull out the candidate pathway given the gene, colocalised
# metabolite, and the metabolite that shares a pathway with colocalised
# metabolite.
#
# gene = gene name
# hmdb.met = colocalised metabolite (HMDB name)
# path.met = metabolite that is part of the pathway shared with the colocalised
#            metabolite
get_cand_path = function(gene, hmdb.met, path.met) {
  met = met_prot %>% filter(gene_name == gene) %>% pull(metabolite) %>% unique
  met = unique(c(met, hmdb.met, path.met))
  all_path = map_dfr(met, ~ met_path %>% filter(metabolite == .x))
  # Need to focus on the pathways that contain both the colocalised metabolite
  # and also the metabolite affected by the gene at the gene's locus
  hmdb.met_path = met_path %>% filter(metabolite == hmdb.met) %>% pull(name) %>% unique
  path.met_path = met_path %>% filter(metabolite == path.met) %>% pull(name) %>% unique
  all_path = all_path %>% filter(name %in% hmdb.met_path, name %in% path.met_path)
  all_path = all_path %>% filter(metabolite %in% c(hmdb.met, path.met))
  all_path = all_path %>% pull(name) %>% unique
  res = data.frame(gene_name = rep(gene, length(all_path)), metabolite.hmdb = rep(hmdb.met, length(all_path)), metabolite.path = rep(path.met, length(all_path)), all_path)
  return(res)
}

# Pull out candidate pathway(s) (and metabolites) of all the genes that
# affected at least one metabolite that is part of the pathway that includes
# the colocalised metabolite
dat_prot_sub = dat_prot %>% distinct(gene_name, metabolite.hmdb, metabolite.pathway)
cand_pathway = pmap_dfr(list(dat_prot_sub$gene_name, dat_prot_sub$metabolite.hmdb, dat_prot_sub$metabolite.pathway), ~ get_cand_path(..1, ..2, ..3))
cand_pathway = cand_pathway %>% distinct

head(dat_prot) %>% as.data.frame
dat_prot_path = left_join(dat_prot, cand_pathway, by = c('metabolite.hmdb', 'metabolite.pathway' = 'metabolite.path', 'gene_name'), relationship = 'many-to-many')

# TODO: from here

dat_prot_path = dat_prot_path %>% select()

tmp = dat_prot_path %>% distinct(metabolite.hmdb, metabolite.pathway, gene_name, all_path) %>% group_by(gene_name, all_path) %>% summarize(total = n()) %>% arrange(desc(total)) %>% slice(1:3)

# tmp = dat_prot_path %>% distinct(metabolite.hmdb, metabolite.pathway, gene_name, all_path) %>% group_by(gene_name, all_path) %>% summarize(total = n()) %>% filter(gene_name == 'ACAD11') %>% arrange(desc(total))
head(tmp)
tmp %>% filter(gene_name == 'DGAT2')
tmp %>% filter(gene_name == 'MOGAT2')


head(dat_prot_path) %>% as.data.frame
dim(dat_prot_path)
head(dat_prot)
dim(dat_prot)
head(cand_pathway)
dim(cand_pathway)

dim(top_pathway)
length(unique(dat_prot$name))
length(unique(top_pathway$name))
head(top_pathway)

head(dat_prot) %>% as.data.frame
dat_prot %>% filter(gene_name == 'CPS1') %>% as.data.frame
dat_prot %>% filter(gene_name == 'ODC1') %>% as.data.frame


test = top_pathway %>% distinct(name, gene_name )
test = left_join(dat_prot, top_pathway, by = c('metabolite.pathway' = 'metabolite'), relationship = 'many-to-many')
head(test) %>% as.data.frame
head(dat_prot)
dim(test)

met_path %>% filter(name == 'Dimethylglycine Dehydrogenase Deficiency')
head(met_path)

test2 = met_path %>% filter(metabolite == 'methionine sulfoxide') %>% pull(name)
test3 = met_path %>% filter(name %in% test2)
test3 = inner_join(test3, met_prot, by = 'metabolite', relationship = 'many-to-many')
head(test3)
test3 %>% filter(gene_name == 'PRIM1') %>% select(metabolite, gene_name, name.x, name.y)
dat_prot %>% filter(metabolite.hmdb == 'methionine sulfoxide', gene_name == 'PRIM1')


head(dat_prot)
test = left_join(dat_prot, top_pathway, by = c('metabolite.pathway' = 'metabolite', 'gene_name'), relationship = 'many-to-many') %>% select(locus:name.protein, gene_name, name) %>% distinct
test %>% filter(gene_name == 'AGMAT')
test %>% filter(gene_name == 'ALDH3A1')
test %>% filter(gene_name == 'WDR1')

met_prot %>% filter(gene_name == 'ALDH3A1') %>% pull(metabolite)
met_path %>% filter(gene_name == 'ALDH3A1')

length(unique(test$name))

dim(dat_prot)
dim(test)
head(test)

top_pathway %>% filter(gene_name == 'GLS2')

write.table(test, 'results/hmdb_match/coloc_results.gene_pathway.txt', sep = '\t', col.names = T, row.names = F, quote = F)

length(unique(test$name))
head(test)
dim(test)

head(top_pathway)
top_pathway %>% distinct %>% dim
top_pathway %>% filter(gene_name == 'RRM2')
top_pathway %>% filter(gene_name == 'SLC17A1') %>% pull(metabolite) %>% unique
top_pathway %>% filter(gene_name == 'SLC17A1') %>% distinct(name, gene_name) %>%  head(., 10)
top_pathway  %>% distinct(gene_name, name) %>%  group_by(gene_name) %>% summarize(count =n()) %>% arrange(desc(count))
dim(top_pathway)
dat_prot %>% filter(metabolite.pathway == 'sodium') %>% pull(gene_name) %>% unique
dat_prot %>% filter(metabolite.pathway == 'phosphate') %>% pull(gene_name) %>% unique



