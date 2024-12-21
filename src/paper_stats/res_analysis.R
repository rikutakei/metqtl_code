library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(VennDiagram)
library(RColorBrewer)

setwd('..')

source('src/paper_stats/functions.R')

# Load all the coloc results
metsim = read.table('results/coloc_res/metsim/merged_res/coloc_merged.pp0.8.txt', sep = '\t', header = T, quote = '')
chen = read.table('results/coloc_res/chen/merged_res/coloc_merged.pp0.8.txt', sep = '\t', header = T, quote = '')
schl_plasma = read.table('results/coloc_res/schlosser/plasma/merged_res/coloc_merged.pp0.8.plasma.txt', sep = '\t', header = T, quote = '')
schl_urine = read.table('results/coloc_res/schlosser/urine/merged_res/coloc_merged.pp0.8.urine.txt', sep = '\t', header = T, quote = '')

# Only focus on full cohort results (276 loci):
metsim = metsim %>% filter(cohort.locus == 'full')
chen = chen %>% filter(cohort.locus == 'full')
schl_plasma = schl_plasma %>% filter(cohort.locus == 'full')
schl_urine = schl_urine %>% filter(cohort.locus == 'full')

metsim %>% group_by(locus) %>% summarise(count = n()) %>% arrange(desc(count))

length(unique(metsim$compound))
length(unique(chen$compound))
length(unique(schl_plasma$compound))
length(unique(schl_urine$compound))

# Add chr/pos info:
loci = read.table('data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% select(SNP:pos) %>% distinct

metsim = left_join(metsim, loci)
chen = left_join(chen, loci)
schl_plasma = left_join(schl_plasma, loci)
schl_urine = left_join(schl_urine, loci)

write.table(metsim, 'results/supp_table/metsim_coloc.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(chen, 'results/supp_table/chen_coloc.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(schl_plasma, 'results/supp_table/schl_plasma_coloc.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(schl_urine, 'results/supp_table/schl_urine_coloc.txt', sep = '\t', col.names = T, row.names = F, quote = F)

################################################################################
# Take a look at per-locus colocalisation results
# - how many metabolites colocalised per locus?
# - average number of metabolites?
# - what locus had the most/least metabolites?

metsim_perloc = per_locus(metsim, metsim = T)
chen_perloc = per_locus(chen)
schl_plasma_perloc = per_locus(schl_plasma)
schl_urine_perloc = per_locus(schl_urine)

table(metsim_perloc$count)
table(chen_perloc$count)
table(schl_plasma_perloc$count)
table(schl_urine_perloc$count)

length(unique(metsim_perloc$locus))
length(unique(chen_perloc$locus))
length(unique(schl_plasma_perloc$locus))
length(unique(schl_urine_perloc$locus))

mean(metsim_perloc$count)
mean(chen_perloc$count)
mean(schl_plasma_perloc$count)
mean(schl_urine_perloc$count)

# What is the breakdown of the categories of all METSIM metabolites:
all_metsim = read.table('data/metqtl_pheno/metsim_phenotypes.tsv', sep = '\t', header = T, stringsAsFactors = F)
all_cat = all_metsim %>% group_by(category) %>% summarize(count = n()) %>% mutate(prop = count / sum(count))
all_cat = add_color(all_cat)

all_cat = all_cat %>% mutate(category = gsub('_', ' ', category), locus = 'All METSIM metabolites (n = 1,391)', broad = 'All METSIM metabolites (n = 1,391)')

all_cat_plot = all_cat %>% ggplot(aes(x = broad, y = prop, fill = category))
all_cat_plot = all_cat_plot +
  scale_fill_manual(values = unique(all_cat$color), labels = sort(unique(all_cat$category)), name = 'Category') +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  labs(x = '', y = 'Proportion of colocalised metabolite categories at each locus') +
  geom_col()

# all_cat_plot

write.table(all_cat, 'results/supp_table/metsim_per_category_summary.total.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# What is the breakdown of the categories of 633 colocalised metabolites in
# METSIM:
cat_sub = metsim %>% distinct(compound, category)
cat_sub = cat_sub %>% group_by(category) %>% summarize(count = n()) %>% mutate(prop = count / sum(count))
cat_sub = add_color(cat_sub)

cat_sub = cat_sub %>% mutate(category = gsub('_', ' ', category), locus = 'Colocalized metabolites (n = 633)', broad = 'colocalized metabolites (n = 633)')

cat_plot = cat_sub %>% ggplot(aes(x = broad, y = prop, fill = category))
cat_plot = cat_plot +
  scale_fill_manual(values = unique(cat_sub$color), labels = sort(unique(cat_sub$category)), name = 'Category') +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  labs(x = '', y = 'Proportion of colocalised metabolite categories at each locus') +
  geom_col()

# cat_plot

write.table(cat_sub, 'results/supp_table/metsim_per_category_summary.633colocalized.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# NOTE: Below analysis may or may not be included in the final version (Nov 20
# 2024)
#
# For the METSIM dataset, plot the proportion of the metabolite categories
# influenced by the loci:
#
# NOTE: Restrict it to those that had >= 5 (or 10) metabolites at the locus, since
# majority of the loci colocalised with 3 or less metabolites

# plot_metsim = metsim_perloc %>% filter(count >= 5) %>% pivot_longer(cols = Amino_Acid:last_col(), names_to = 'category', values_to = 'met_count') %>% mutate(prop = met_count / count)
plot_metsim = metsim_perloc %>% filter(count >= 10) %>% pivot_longer(cols = Amino_Acid:last_col(), names_to = 'category', values_to = 'met_count') %>% mutate(prop = met_count / count)

per_cat_summary = plot_metsim %>% group_by(category) %>% summarize(count = sum(met_count), mean_prop = mean(prop), max_prop = max(prop), min_prop = min(prop), sd_prop = sd(prop)) %>% arrange(desc(mean_prop))
per_cat_summary
write.table(per_cat_summary, 'results/supp_table/metsim_per_category_summary.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# plot_metsim %>% filter(SNP == 'rs2231142', prop > 0) %>% as.data.frame
# plot_metsim %>% filter(prop > 0, category == 'Xenobiotics') %>% as.data.frame

plot_metsim = add_color(plot_metsim)

# Rename the locus for plotting:
plot_metsim = plot_metsim %>% separate(locus, sep = '_', into = c('chr', 'start', 'end'), remove = F) %>% unite(col = 'range', start, end, sep = '-', remove = F) %>% arrange(str_sort(locus, numeric = T))

loci_gene = read.table('src/paper_stats/loci_genes.txt', sep = '\t', header = T, stringsAsFactors = F)
loci_gene = loci_gene %>% select(-locus.plot)

plot_metsim = left_join(plot_metsim, loci_gene, by = 'SNP')
plot_metsim = plot_metsim %>% unite(col = 'locus.plot', chr, range, sep = ':', remove = F) %>% mutate(locus.plot = paste(gsub(' /', ';', SNP), ' / ', genes, ' (', locus.plot, ')', sep = '')) %>% mutate(locus.plot = gsub('/ NA ', '', locus.plot))

plot_level = plot_metsim %>% distinct(locus, locus.plot) %>% as.data.frame
rownames(plot_level) = plot_level$locus
plot_level = plot_level[str_sort(plot_level$locus, numeric = T),]
plot_level = plot_level$locus.plot

plot_metsim$locus.plot = factor(plot_metsim$locus.plot, levels = plot_level)
plot_metsim$locus = plot_metsim$locus.plot

plot_perloc(plot_metsim, facet = F, outname = 'results/plots/coloc/metsim_category_prop.stacked.png')
plot_perloc(plot_metsim, facet = T, outname = 'results/plots/coloc/metsim_category_prop.facet.png')

# Combine the overal category data with the per-locus data:
overall = plot_metsim %>% mutate(broad = 'Per-locus') %>% select(colnames(cat_sub))
overall = rbind(overall, cat_sub)
overall = rbind(overall, all_cat)
overall = overall %>% mutate(category = gsub('_', ' ', category))

overall_plot = overall %>% ggplot(aes(x = locus, y = prop, fill = category)) +
  facet_grid( ~ broad, scales = 'free', space = 'free') +
  scale_fill_manual(values = unique(overall$color), labels = sort(unique(overall$category)), name = 'Category') +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x = '', y = 'Proportion of colocalised metabolite categories at each locus') +
  geom_col()

overall_plot

ggsave('results/plots/coloc/metsim_category_prop.stacked.total.png', width = 10, height = 8)

# GLS2 locus also seems suspicious in terms of metabolite category ratio:
# plot_metsim %>% filter(SNP == "rs2638315") %>% as.data.frame

##############################################################################
# TODO: for the proportion of the categories, is this what you would expect to
# see based on the original data's proportion for each category? e.g. if the
# METSIM data metabolites are 70% lipids, then a locus with 75% lipids may not
# as interesting as it may seem
##############################################################################

# Repeat this process with the Schlosser data set (for those that are present
# in the METSIM data set)
metsim_hmdb = read.table('results/hmdb_match/metsim/all_metsim.hmdb_merge.txt', sep = '\t', quote = '', header = T, stringsAsFactors = F)
metsim_hmdb = metsim_hmdb %>% select(category, metabolite.hmdb)

# Merge the HMDB metabolite names to Schlosser data set:
schl_hmdb = read.table('results/hmdb_match/schlosser/coloc_results.hmdb_merge.match_only.txt', sep = '\t', header = T, stringsAsFactors = F)
schl_hmdb = schl_hmdb %>% select(compound, metabolite.hmdb, category)

schl_plasma_hmdb = left_join(schl_plasma, schl_hmdb, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb))
schl_urine_hmdb = left_join(schl_urine, schl_hmdb, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb))

# Add category information to Schlosser based on the METSIM-HMDB metabolite matches
schl_plasma_hmdb = left_join(schl_plasma_hmdb, metsim_hmdb, by = 'metabolite.hmdb', relationship = 'many-to-many', suffix = c('.schl', '')) %>% select(locus:PP.H4.abf, category)
schl_urine_hmdb = left_join(schl_urine_hmdb, metsim_hmdb, by = 'metabolite.hmdb', relationship = 'many-to-many', suffix = c('.schl', '')) %>% select(locus:PP.H4.abf, category)

# Remove those that didn't match:
schl_plasma_hmdb = schl_plasma_hmdb %>% filter(!is.na(category)) %>% distinct
schl_urine_hmdb = schl_urine_hmdb %>% filter(!is.na(category)) %>% distinct

# Group the categories per locus
schl_plasma_hmdb = per_locus(schl_plasma_hmdb, metsim = T)
schl_urine_hmdb = per_locus(schl_urine_hmdb, metsim = T)

# Now calculate the proportion of each category and plot it as in METSIM
plot_schl_plasma = schl_plasma_hmdb %>% filter(count >= 10) %>% pivot_longer(cols = Amino_Acid:last_col(), names_to = 'category', values_to = 'met_count') %>% mutate(prop = met_count / count)
plot_schl_urine = schl_urine_hmdb %>% filter(count >= 10) %>% pivot_longer(cols = Amino_Acid:last_col(), names_to = 'category', values_to = 'met_count') %>% mutate(prop = met_count / count)

plot_schl_plasma = add_color(plot_schl_plasma)
plot_schl_urine = add_color(plot_schl_urine)

plot_schl_plasma %>% filter(SNP == 'rs2231142') %>% pull(metabolites) %>% unique %>% strsplit(., '; ')

plot_perloc(plot_schl_plasma, facet = F, outname = 'results/plots/coloc/schl_plasma_category_prop.stacked.png')
plot_perloc(plot_schl_plasma, facet = T, outname = 'results/plots/coloc/schl_plasma_category_prop.facet.png')

plot_perloc(plot_schl_urine, facet = F, outname = 'results/plots/coloc/schl_urine_category_prop.stacked.png')
plot_perloc(plot_schl_urine, facet = T, outname = 'results/plots/coloc/schl_urine_category_prop.facet.png')

################################################################################
# Some stats for the paper:
# - number of Schlosser metabolites that were able to be matched with METSIM, based on HMDB names
length(unique(schl_plasma$compound))
length(unique(schl_urine$compound))
length(unique(unlist(map(schl_plasma_hmdb$compound_id, ~ strsplit(.x, '; ')))))
length(unique(unlist(map(schl_urine_hmdb$compound_id, ~ strsplit(.x, '; ')))))

# Take a look at the overlap of METSIM, Chen, and Schlosser plasma:
metsim_hmdb = read.table('results/hmdb_match/metsim/all_metsim.hmdb_merge.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
metsim_hmdb = metsim_hmdb %>% select(phenostring:match.type)
metsim_hmdb = left_join(metsim, metsim_hmdb, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb)) %>% mutate(metabolite = gsub('"', '', metabolite), metabolite.hmdb = gsub('"', '', metabolite.hmdb), synonym = gsub('"', '', synonym))

chen_hmdb = read.table('results/hmdb_match/chen/coloc_results.hmdb_merge.txt', sep = '\t', header = T, stringsAsFactors = F)
chen_hmdb = chen_hmdb %>% select(phenostring:match.type)
chen_hmdb = left_join(chen, chen_hmdb, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb)) %>% mutate(metabolite = gsub('"', '', metabolite), metabolite.hmdb = gsub('"', '', metabolite.hmdb), synonym = gsub('"', '', synonym))

schl_hmdb = read.table('results/hmdb_match/schlosser/coloc_results.hmdb_merge.match_only.txt', sep = '\t', header = T, stringsAsFactors = F)
schl_hmdb = schl_hmdb %>% select(compound, category:match.type)

schl_plasma_hmdb = left_join(schl_plasma, schl_hmdb, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb)) %>% mutate(metabolite = gsub('"', '', metabolite), metabolite.hmdb = gsub('"', '', metabolite.hmdb), synonym = gsub('"', '', synonym))
schl_urine_hmdb = left_join(schl_urine, schl_hmdb, relationship = 'many-to-many') %>% filter(!is.na(metabolite.hmdb)) %>% mutate(metabolite = gsub('"', '', metabolite), metabolite.hmdb = gsub('"', '', metabolite.hmdb), synonym = gsub('"', '', synonym))

# Count up some numbers
length(unique(metsim$phenostring))
length(unique(chen$phenostring))
length(unique(schl_plasma$phenostring))
length(unique(schl_urine$phenostring))

# NOTE: the numbers for *_hmdb are slightly 'wrong' when you compare between
# data sets, since phenostring may be different between data sets (e.g. one
# data set may have two unique phenostring for one metabolite.hmdb, whereas
# another may only have one)
length(unique(metsim_hmdb$phenostring))
length(unique(chen_hmdb$phenostring))
length(unique(schl_plasma_hmdb$phenostring))
length(unique(schl_urine_hmdb$phenostring))

# Since there may be more than one HMDB metabolite associated with each
# phenotype, it is difficult to make a Venn diagram with the phenostring names
# from different data sets.
#
# So start by merging all three data sets together based on the HMDB
# metabolites, then assign one representative (cleaned up) phenostring for each
# metabolite:
merged = metsim_hmdb %>% select(metabolite, metabolite.hmdb)
merged = chen_hmdb %>% select(metabolite, metabolite.hmdb) %>% full_join(merged, ., by = 'metabolite.hmdb', suffix = c('', '.chen'), relationship = 'many-to-many')
merged = schl_urine_hmdb %>% select(metabolite, metabolite.hmdb) %>% full_join(merged, ., by = 'metabolite.hmdb', suffix = c('', '.schl_urine'), relationship = 'many-to-many')
merged = schl_plasma_hmdb %>% select(metabolite, metabolite.hmdb) %>% full_join(merged, ., by = 'metabolite.hmdb', suffix = c('.metsim', '.schl_plasma'), relationship = 'many-to-many') %>% distinct

master_phenostring = merged %>% pivot_longer(cols = c('metabolite.metsim', 'metabolite.chen', 'metabolite.schl_urine', 'metabolite.schl_plasma'), names_to = 'study', values_to = 'metabolite') %>% drop_na %>% group_by(metabolite.hmdb) %>% slice(1) %>% ungroup %>% distinct(metabolite.hmdb, metabolite)

head(master_phenostring)
tmp = master_phenostring %>% filter(metabolite == "undecenoylcarnitine") %>% pull(metabolite.hmdb)
head(metsim_hmdb)
metsim_hmdb %>% filter(metabolite.hmdb %in% tmp )

write.table(master_phenostring, 'results/hmdb_match/master_phenostring.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Using the master phenostring data, add back the phenostring to the data sets,
# then look at the overlap (which should now be the same across data sets):
overlap_metsim = left_join(metsim_hmdb, master_phenostring, by = 'metabolite.hmdb', suffix = c('', '.common'))
overlap_chen = left_join(chen_hmdb, master_phenostring, by = 'metabolite.hmdb', suffix = c('', '.common'))
overlap_schl = left_join(schl_plasma_hmdb, master_phenostring, by = 'metabolite.hmdb', suffix = c('', '.common'))

# Number of metabolites that had HMDB matches:
length(unique(overlap_metsim$metabolite.common))
length(unique(overlap_chen$metabolite.common))
length(unique(overlap_schl$metabolite.common))

overlap_list = list(unique(overlap_metsim$metabolite.common), unique(overlap_chen$metabolite.common), unique(overlap_schl$metabolite.common))

venn.diagram(x = overlap_list,
             category.names = c('METSIM', 'CLSA', 'GCKD'),
             filename = 'results/plots/venn/plasma_metabolite_overlap_venn.png',
             output = T,
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = brewer.pal(3, 'Set1'),
             # Numbers
             cex = 1.6,
             fontfamily = "sans",
             # Set names
             cat.cex = 1.2,
             cat.fontface = "bold",
             # cat.default.pos = "outer",
             # cat.pos = c(-27, 27, 135),
             # cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1
)

venn.diagram(x = overlap_list,
             category.names = c('METSIM', 'CLSA', 'GCKD'),
             filename = 'results/plots/venn/plasma_metabolite_overlap_venn.bw.png',
             output = T,
             # Numbers
             cex = 1.6,
             fontfamily = "sans",
             # Set names
             cat.cex = 1.2,
             cat.fontface = "bold",
             # cat.default.pos = "outer",
             # cat.pos = c(-27, 27, 135),
             # cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1
)

# Save the list of 141 metabolites common in all three data sets:
all_common = table(unlist(overlap_list))
all_common = names(all_common)[which(all_common == 3)]

tmp_metsim = overlap_metsim %>% select(metabolite.common, SNP, locus)
tmp_chen = overlap_chen %>% select(metabolite.common, SNP, locus)
tmp_schl = overlap_schl %>% select(metabolite.common, SNP, locus)
tmp_all_common = rbind(tmp_metsim, tmp_chen)
tmp_all_common = rbind(tmp_all_common, tmp_schl)

tmp_all_common = tmp_all_common %>% filter(metabolite.common %in% all_common) %>% group_by(metabolite.common) %>% summarise(coloc_loci = paste(unique(locus), collapse = '; '), coloc_snp = paste(unique(SNP), collapse = '; '))

write.table(tmp_all_common, 'results/met_overlap/all_plasma_overlap.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Look at the overlap between 141 plasma metabolites in urine
#
# First check how many urine metabolites matched with HMDB (merge with master
# phenostring, then count unique metabolites):
overlap_schl_urine = left_join(schl_urine_hmdb, master_phenostring, by = 'metabolite.hmdb', suffix = c('', '.common'))

# 218 urine metabolites matched with HMDB
length(unique(overlap_schl_urine$metabolite.common))

# Now see how many overlapped with the 141 plasma metabolites
# (Out of 218 urine metabolites, 52 overlapped with plasma metabolites)
plasma_urine_overlap = overlap_schl_urine %>% filter(metabolite.common %in% all_common) %>% pull(metabolite.common) %>% unique
length(plasma_urine_overlap)

tmp_plasma_urine_overlap = tmp_all_common %>% filter(metabolite.common %in% plasma_urine_overlap)

# Save the list of 52 common metabolites between plasma and urine:
write.table(tmp_plasma_urine_overlap, 'results/met_overlap/plasma_urine_overlap.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Save the list of 89 metabolites unique in plasma:
unique_plasma = all_common[!(all_common %in% plasma_urine_overlap)]
length(unique_plasma)

tmp_unique_plasma = tmp_all_common %>% filter(metabolite.common %in% unique_plasma)

write.table(tmp_unique_plasma, 'results/met_overlap/unique_plasma.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(tmp_unique_plasma, 'results/supp_table/unique_plasma.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Generate a Venn diagram of the 141 plasma vs 218 urine metabolites:
urine_plasma_overlap = list(unique(all_common), unique(overlap_schl_urine$metabolite.common))

venn.diagram(x = urine_plasma_overlap,
             category.names = c('Plasma', 'Urine'),
             filename = 'results/plots/venn/plasma_urine_metabolite_overlap_venn.png',
             output = T,
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = brewer.pal(3, 'Set1')[1:2],
             # Numbers
             cex = 1.6,
             fontfamily = "sans",
             # Set names
             cat.cex = 1.2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-180, -180),
             cat.dist = c(0.025, 0.025),
             cat.fontfamily = "sans"
)

venn.diagram(x = urine_plasma_overlap,
             category.names = c('Plasma', 'Urine'),
             filename = 'results/plots/venn/plasma_urine_metabolite_overlap_venn.bw.png',
             output = T,
             # Numbers
             cex = 1.6,
             fontfamily = "sans",
             # Set names
             cat.cex = 1.2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-180, -180),
             cat.dist = c(0.025, 0.025),
             cat.fontfamily = "sans"
)

################################################################################
# List up all the loci/SNPs involved with the 141 metabolites for pathway
# enrichment analysis
metsim_snp = overlap_metsim %>% filter(metabolite.common %in% all_common) %>% pull(SNP) %>% unique
chen_snp = overlap_chen %>% filter(metabolite.common %in% all_common) %>% pull(SNP) %>% unique
schl_snp = overlap_schl %>% filter(metabolite.common %in% all_common) %>% pull(SNP) %>% unique

all_snp = unique(c(metsim_snp, chen_snp, schl_snp))

length(all_snp)

writeLines(all_snp, 'results/met_overlap/common_plasma_metabolite_snps.txt')

# List up those that were present only in plasma
plasma_met = all_common[which(!(all_common %in% plasma_urine_overlap))]
metsim_snp = overlap_metsim %>% filter(metabolite.common %in% plasma_met) %>% pull(SNP) %>% unique
chen_snp = overlap_chen %>% filter(metabolite.common %in% plasma_met) %>% pull(SNP) %>% unique
schl_snp = overlap_schl %>% filter(metabolite.common %in% plasma_met) %>% pull(SNP) %>% unique

plasma_only = unique(c(metsim_snp, chen_snp, schl_snp))

length(plasma_only)

writeLines(plasma_only, 'results/met_overlap/plasma_only_snps.txt')

# And those that overlapped with urine
metsim_snp = overlap_metsim %>% filter(metabolite.common %in% plasma_urine_overlap) %>% pull(SNP) %>% unique
chen_snp = overlap_chen %>% filter(metabolite.common %in% plasma_urine_overlap) %>% pull(SNP) %>% unique
schl_snp = overlap_schl %>% filter(metabolite.common %in% plasma_urine_overlap) %>% pull(SNP) %>% unique
schl_urine_snp = overlap_schl %>% filter(metabolite.common %in% plasma_urine_overlap) %>% pull(SNP) %>% unique

plasma_urine_snp = unique(c(metsim_snp, chen_snp, schl_snp, schl_urine_snp))

length(plasma_urine_snp)

writeLines(plasma_urine_snp, 'results/met_overlap/plasma_urine_snps.txt')

################################################################################
# Summarise the stats for the loci that colocalised per compound/metabolite:
metsim_percomp = per_compound(metsim)
chen_percomp = per_compound(chen)
schl_plasma_percomp = per_compound(schl_plasma)
schl_urine_percomp = per_compound(schl_urine)

table(metsim_percomp$count)
table(chen_percomp$count)
table(schl_plasma_percomp$count)
table(schl_urine_percomp$count)

length(unique(metsim_percomp$phenostring))
length(unique(chen_percomp$phenostring))
length(unique(schl_plasma_percomp$phenostring))
length(unique(schl_urine_percomp$phenostring))

mean(metsim_percomp$count)
mean(chen_percomp$count)
mean(schl_plasma_percomp$count)
mean(schl_urine_percomp$count)

# Again, only for the METSIM data, take a look at the number of loci
# colocalised based on the category of metabolites
#
# NOTE: the plot doesn't tell you a lot - just that most metabolites are
# affected/colocalise with 1 locus

plot_metcomp = metsim %>% distinct(compound, category) %>% left_join(metsim_percomp, .)
plot_metcomp = add_color(plot_metcomp)

plot_metcomp %>% ggplot(aes(x = count, fill = category)) +
  scale_fill_manual(values = unique(plot_metcomp$color), labels = sort(unique(plot_metcomp$category)), name = 'Category') +
  facet_wrap( ~ category) +
  theme_bw() +
  geom_bar()

metsim_metcomp = metsim %>% distinct(compound, category) %>% left_join(metsim_percomp, .) %>% select(compound, phenostring, category, count:loci)

################################################################################

# OLD STUFF (I think...):
# How many of the metabolites from METSIM also colocalised in Schlosser plasma
# data set?

# First need to make the metabolite names (somewhat) consistent, so use the
# HMDB-matched data set

metsim_hmdb = read.table("results/hmdb_match/metsim/coloc_results.hmdb_merge.txt", sep = '\t', header = T, stringsAsFactors = F)
schl_hmdb = read.table("results/hmdb_match/schlosser/coloc_results.hmdb_merge.txt", sep = '\t', header = T, stringsAsFactors = F) %>% filter(category == 'plasma')

# Number of metabolites that were able to match with HMDB metabolites:
# METSIM: 390/486 (80.2%)
# Schlosser plasma: 307/480 (64.0%)
#
# NOTE: these are the numbers based on the "cleaned" metabolite names, so the
# lipids will have been simplified

length(unique(metsim_hmdb$metabolite))
metsim_met = metsim_hmdb %>% filter(!is.na(match.type)) %>% pull(metabolite.hmdb) %>% unique
length(metsim_met)

length(unique(schl_hmdb$metabolite))
schl_met = schl_hmdb %>% filter(!is.na(match.type)) %>% pull(metabolite.hmdb) %>% unique
length(schl_met)

# How many are in common? (193 in common out of a total of 504)
#
# So, there are 197 specific to METSIM and 114 specific to Schlosser plasma
common = metsim_met[which(metsim_met %in% schl_met)]
total = unique(c(metsim_met, schl_met))

length(common)
length(total)

writeLines(common, 'results/met_overlap/metsim_schlosser_plasma_overlap.txt')

# Let's see what category these metabolites are in (based on the METSIM
# category data)
common_met = metsim_hmdb %>% filter(metabolite.hmdb %in% common) %>% select(locus, compound:SNP, metabolite:match.type)

common_met %>% distinct(category, phenostring) %>% mutate(category = gsub('_', ' ', category)) %>% ggplot(aes(x = category, fill = category)) +
  scale_fill_manual(values = color$color, name = 'Category') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = 'Categories of metabolites that colocalised with gout in both METSIM and Schlosser (plasma) data sets', x = '', y = 'Number of metabolites') +
  geom_bar(show.legend = F)

ggsave('results/plots/coloc/common_met_category.png', width = 10, height = 11)

################################################################################

# OLD STUFF:
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

