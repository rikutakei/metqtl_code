library(dplyr)
library(tidyr)

setwd('..')

# Load gene UCSC information
gene = read.table('src/locuszoom/LocusZooms/UCSC_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactor = F)
gene = gene %>% filter(!(Chrom %in% c('chrM', 'chrY'))) %>% mutate(Chrom =gsub('X', '23', Chrom), Chrom = as.integer(gsub('chr', '', Chrom)))

# Make data for the SNPs involved with the metabolites
loci = read.table('data/gwas/indep_snps_from_paper.eur_only.txt', sep = '\t', header = T, stringsAsFactor = F)
loci = loci %>% select(SNP, chr, pos) %>% distinct

all_common = readLines('results/met_overlap/common_plasma_metabolite_snps.txt')
plasma_only = readLines('results/met_overlap/plasma_only_snps.txt')
plasma_urine = readLines('results/met_overlap/plasma_urine_snps.txt')

all_common = loci %>% filter(SNP %in% all_common)
plasma_only = loci %>% filter(SNP %in% plasma_only)
plasma_urine = loci %>% filter(SNP %in% plasma_urine)

all_common = all_common %>% mutate(start = pos - 500000, end = pos + 500000)
plasma_only = plasma_only %>% mutate(start = pos - 500000, end = pos + 500000)
plasma_urine = plasma_urine %>% mutate(start = pos - 500000, end = pos + 500000)

# Overlap the gene information with the SNP +/-500kb regions:
all_common = gene %>% select(Chrom, Start, End, Gene) %>% left_join(all_common, ., by = c('chr' = 'Chrom'), relationship = 'many-to-many')
all_common = all_common %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End))
all_common_genes = unique(all_common$Gene)

plasma_only = gene %>% select(Chrom, Start, End, Gene) %>% left_join(plasma_only, ., by = c('chr' = 'Chrom'), relationship = 'many-to-many')
plasma_only = plasma_only %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End))
plasma_only_genes = unique(plasma_only$Gene)

plasma_urine = gene %>% select(Chrom, Start, End, Gene) %>% left_join(plasma_urine, ., by = c('chr' = 'Chrom'), relationship = 'many-to-many')
plasma_urine = plasma_urine %>% filter(between(Start, start, end) | between(End, start, end) | between(start, Start, End))
plasma_urine_genes = unique(plasma_urine$Gene)

writeLines(all_common_genes, 'results/pathway/all_common_met_genes.all.txt')
writeLines(plasma_only_genes, 'results/pathway/plasma_only_met_genes.all.txt')
writeLines(plasma_urine_genes, 'results/pathway/plasma_urine_met_genes.all.txt')

# Note that these are all genes, including pseudogenes, ncRNA, etc., so make
# a separate list of just protein coding genes:
prot_only = gene %>% filter(Coding == 'proteincoding') %>% pull(Gene)

all_common_genes_prot = all_common_genes[which(all_common_genes %in% prot_only)]
plasma_only_genes_prot = plasma_only_genes[which(plasma_only_genes %in% prot_only)]
plasma_urine_genes_prot = plasma_urine_genes[which(plasma_urine_genes %in% prot_only)]

writeLines(all_common_genes_prot, 'results/pathway/all_common_met_genes.prot_only.txt')
writeLines(plasma_only_genes_prot, 'results/pathway/plasma_only_met_genes.prot_only.txt')
writeLines(plasma_urine_genes_prot, 'results/pathway/plasma_urine_met_genes.prot_only.txt')

