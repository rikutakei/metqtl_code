# Analyze colocalised metabolites that have matching HMDB metabolites
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)

setwd('../..')

# Load in list of 144 common plasma metabolites
met = readLines('results/met_overlap/all_plasma_overlap.txt')

# Pull out the HMDB names from the coloc results.
#
# NOTE: Remember that these 144 metabolites are the cleaned up names of the
# phenotypes from all three data sets (METSIM, Chen, and Schlosser), and
# multiple HMDB metabolite names can match to a single phenostring
metsim = read.table('results/hmdb_match/metsim/coloc_results.hmdb_merge.match_only.txt', sep = '\t', header = T, stringsAsFactors = F)
schl = read.table('results/hmdb_match/schlosser/coloc_results.hmdb_merge.match_only.txt', sep = '\t', header = T, stringsAsFactors = F)
chen = read.table('results/hmdb_match/chen/coloc_results.hmdb_merge.match_only.txt', sep = '\t', header = T, stringsAsFactors = F)

dat = map_dfr(list(metsim, schl, chen), ~ .x %>% select(SNP, phenostring, metabolite:match.type) %>% distinct)
dat = dat %>% distinct

# Merge the master phenostring mapping data set
master = read.table('results/hmdb_match/master_phenostring.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

dat = left_join(dat, master, by = 'metabolite.hmdb', suffix = c('.original', ''))

# Pull out the target metabolites for lookup
# (160 HMDB metabolites to look up pathways for)
target = dat %>% filter(metabolite %in% met)
target = target %>% pull(metabolite.hmdb) %>% unique

# TODO: add which gene(s) showed up directly in the colocalised metabolite (e.g. GLS2 with vitamin A (I think))
#
# Majority of the colocalised metabolites are NOT directly
# controlled/transported by the protein product of the gene at the locus.
# It is likely that the metabolite is either a stable intermediate or an
# end-product of a pathway, where the gene(s) at the locus affects the
# metabolites before and/or after the colocalised metabolite. Also, HMDB is
# a little "patchy" in terms of transporters (e.g. doesn't have ABCG2 or
# SLC22A9 listed for uric acid)

# See if any of the proteins/genes involved in any of the metabolites that is
# part of any of the pathways with the colocalised metabolites is within the
# locus.

# Load in metabolite/pathway data:
met_path = read.table('data/hmdb/hmdb_biopath.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_path = met_path %>% mutate(metabolite = tolower(metabolite))

# Pull out all the pathways that have one of the metabolites in the list
# (2106 unique pathways involved)
path_list = met_path %>% filter(metabolite %in% target) %>% pull(name) %>% unique

# Now pull out all of the metabolites that are involved in any of the pathways
# listed above
full_list = met_path %>% filter(name %in% path_list)

# 6289 metabolites from 2106 pathways
length(unique(full_list$metabolite))

# Format the table so that each metabolite has a pathway entry, and then a list
# of metabolites involved in the listed pathway:
full_list = full_list %>% select(metabolite, name) %>% left_join(full_list, ., by = 'name', suffix = c('', '.pathway'), relationship = 'many-to-many') %>% distinct

# Exclude metabolites that are "inroganic" (e.g. water):
met_class = read.table('data/hmdb/hmdb_chemtax.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_class = met_class %>% mutate(metabolite = tolower(metabolite))

inorg = met_class %>% filter(kingdom == 'Inorganic compounds') %>% select(metabolite, kingdom)
full_list = left_join(full_list, inorg, by = 'metabolite')
full_list = left_join(full_list, inorg, by = c('metabolite.pathway' = 'metabolite'), suffix = c('.met', '.pth'))
full_list = full_list %>% filter(is.na(kingdom.met) & is.na(kingdom.pth)) %>% select(!contains('kingdom'))

# Remove metabolites that are present in the pathways of >= X% of metabolites.
#
# For example, ATP ("adenosine triphosphate") is a part of at least one pathway
# in 97.8% of the metabolites - the odds of linking ATP with a given metabolite
# is high, so remove them to reduce noise in the end result
perc_met = full_list %>% group_by(metabolite) %>% summarize(total = length(unique(metabolite.pathway))) %>% ungroup %>% mutate(perc_met = total / nrow(.))

full_list = left_join(full_list, perc_met, by = 'metabolite')

# There is a very sharp drop-off in the percentage after 9.3%, where the next
# percentage after that is 60.1%:
perc_met %>% arrange(desc(perc_met)) %>% filter(perc_met >= 0.03) %>% distinct(metabolite, perc_met) %>% as.data.frame
full_list %>% ggplot(aes(perc_met)) +
  geom_density()

# Remove metabolites that are part of >10% of other metabolites' pathway:
hi_perc = full_list %>% filter(perc_met >= 0.1) %>% pull(metabolite) %>% unique
full_list = full_list %>% filter(!(metabolite %in% hi_perc) & !(metabolite.pathway %in% hi_perc))

# Load metabolite-gene data:
met_prot = read.table('data/hmdb/hmdb_metprot.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
met_prot = met_prot %>% mutate(metabolite = tolower(metabolite))

# Now link up the enzymes/proteins/genes involved with each metabolite:
full_list = left_join(full_list, met_prot, by = c('metabolite.pathway' = 'metabolite'), suffix = c('.pathway', '.protein'), relationship = 'many-to-many')

# Make a smaller table for merging back the pathways based on metabolites and
# proteins:
path_prot = full_list %>% select(metabolite, name.pathway, metabolite.pathway, uniprot_id, gene_name) %>% distinct

# Need to make the `full_list` smaller in order to merge it with the coloc
# result (to determine which genes/proteins are within the locus) - this is
# done by ignoring the pathway column in which the metabolite is associated
# with (hence `path_prot` was generated beforehand to merge the pathway
# information back in, when/if it is needed later):
full_list = full_list %>% filter(!is.na(gene_name)) %>% select(metabolite, metabolite.pathway:protein_type) %>% distinct

# 3202 genes identified from 6289 metabolites
length(unique(full_list$gene_name))

# Merge the full_list with the list of colocalised loci that were involved with
# the 144 common plasma metabolites, but first add in chr/pos info
loci = read.table('data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', header = T, stringsAsFactor = F)
loci = loci %>% select(locus, SNP, chr, pos) %>% distinct

dat = left_join(dat, loci, relationship = 'many-to-many')
dat = dat %>% mutate(start = pos - 500000, end = pos + 500000)

dat_prot = left_join(dat, full_list, by = c('metabolite.hmdb' = 'metabolite'), relationship = 'many-to-many') %>% filter(!is.na(gene_name))
dat_prot = dat_prot %>% distinct(locus, SNP, chr, pos, start, end, phenostring, metabolite.hmdb, metabolite.pathway, name.protein, uniprot_id, gene_name, protein_type)

# Add gene position information so I can determine which genes are present at
# the locus
gene = read.table('src/locuszoom/LocusZooms/UCSC_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactor = F)
gene = gene %>% filter(!(Chrom %in% c('chrM', 'chrY'))) %>% mutate(Chrom =gsub('X', '23', Chrom), Chrom = as.integer(gsub('chr', '', Chrom)))

dat_prot = left_join(dat_prot, gene, by = c('gene_name' = 'Gene'), relationship = 'many-to-many')

# Now filter out genes that are not within (or overlap with) the coloc region.
# First, those not on the same chromosome
dat_prot = dat_prot %>% filter(chr == Chrom)

# Check positioning of the gene.
# If the SNP position, locus start, or locus end is within the gene region,
# then the gene (or at least part of the gene) is within the locus
dat_prot = dat_prot %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End)) %>% as.data.frame

# Clean up
dat_prot = dat_prot %>% select(locus:End, Coding)

# How many genes at how many loci?
length(unique(dat_prot$gene_name)) # 179 genes
length(unique(dat_prot$SNP)) # 78 loci

writeLines(unique(dat_prot$gene_name), 'results/pathway/hmdb_genes.txt')
# write.table(dat_prot, 'results/hmdb_match/coloc_results.met_pathway_genes.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# How many genes directly affected any one of the 144 plasma metabolite?
direct_gene = dat_prot %>% filter(metabolite.hmdb == metabolite.pathway) %>% pull(gene_name) %>% unique
direct_snp = dat_prot %>% filter(metabolite.hmdb == metabolite.pathway) %>% pull(SNP) %>% unique

length(direct_gene)
length(direct_snp)
