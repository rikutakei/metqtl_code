library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

setwd('..')

all_met = read.table('results/pathway/david_144met_genes.txt', sep = '\t', header = T, stringsAsFactors = F)
all_met_coding = read.table('results/pathway/david_144met_genes.coding_only.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

plasma_only = read.table('results/pathway/david_plasma_only_genes.txt', sep = '\t', header = T, stringsAsFactors = F)
plasma_only_coding = read.table('results/pathway/david_plasma_only_genes.coding_only.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')

plasma_urine = read.table('results/pathway/david_plasma_urine_genes.txt', sep = '\t', header = T, stringsAsFactors = F)
plasma_urine_coding = read.table('results/pathway/david_plasma_urine_genes.coding_only.txt', sep = '\t', header = T, stringsAsFactors = F)

hmdb = read.table('results/pathway/david_hmdb_genes.txt', sep = '\t', header = T, stringsAsFactors = F)

clean_david = function(dat) {
  colnames(dat) = c("category", "term", "count", "percentage", "pvalue", "genes", "list.total", "pop.hits", "pop.total", "fold.enrichment", "bonferroni", "benjamini", "fdr")
  res = dat %>% select(category:term, genes, fold.enrichment, pvalue, bonferroni:fdr)
  # res = res %>% mutate(term = gsub('.*~', '', term), term = gsub('hsa[^:]*:', '', term), term = str_to_title(term))
  res = res %>% mutate(term = gsub('.*~', '', term), term = gsub('hsa[^:]*:', '', term), term = str_to_sentence(term))
  res = res %>% arrange(fdr)
  return(res)
}

all_met = clean_david(all_met)
all_met_coding = clean_david(all_met_coding)
plasma_only = clean_david(plasma_only)
plasma_only_coding = clean_david(plasma_only_coding)
plasma_urine = clean_david(plasma_urine)
plasma_urine_coding = clean_david(plasma_urine_coding)
hmdb = clean_david(hmdb)

all_met_sig = all_met %>% filter(fdr <= 0.05)
all_met_coding_sig = all_met_coding %>% filter(fdr <= 0.05)
plasma_only_sig = plasma_only %>% filter(fdr <= 0.05)
plasma_only_coding_sig = plasma_only_coding %>% filter(fdr <= 0.05)
plasma_urine_sig = plasma_urine %>% filter(fdr <= 0.05)
plasma_urine_coding_sig = plasma_urine_coding %>% filter(fdr <= 0.05)
hmdb_sig = hmdb %>% filter(fdr <= 0.05)

write.table(all_met, 'results/pathway/david_144met_genes.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(all_met_coding, 'results/pathway/david_144met_genes.coding_only.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(plasma_only, 'results/pathway/david_plasma_only_genes.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(plasma_only_coding, 'results/pathway/david_plasma_only_genes.coding_only.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(plasma_urine, 'results/pathway/david_plasma_urine_genes.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(plasma_urine_coding, 'results/pathway/david_plasma_urine_genes.coding_only.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(hmdb, 'results/pathway/david_hmdb_genes.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# What new pathways were found in the HMDB-based gene set compared to the other
# gene sets:
found_path = unique(c(all_met_coding_sig$term, plasma_only_coding_sig$term, plasma_urine_coding_sig$term))

hmdb_new = hmdb_sig %>% filter(!(term %in% found_path))

write.table(hmdb_new, 'results/pathway/david_hmdb_genes.new.txt', sep = '\t', col.names = T, row.names = F, quote = F)

