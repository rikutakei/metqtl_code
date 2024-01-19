# Helper script to generate LZ, given a compound code, an rsID, and cohort
#
# Rscript rs2231142 C98 full
library(dplyr)
library(tidyr)

source('src/locuszoom/LocusZooms/functions/locus_zoom.R')

genes = read.delim('src/locuszoom/LocusZooms/Gencode_GRCh37_Genes_UniqueList2021.txt', header = T, stringsAsFactors = F)

args = commandArgs(trailingOnly = T)
snp = args[1]
compound = args[2]
cohort = args[3]

gwas_file = gsub('SEX', cohort, 'results/gwas/SEX/SNP.gwas.txt')
gwas_file = gsub('SNP', snp, gwas_file)
gwas = read.table(gwas_file, sep = '\t', header = T, stringsAsFactors = F)

met_file = gsub('COM', compound, 'results/metqtl/COM/SNP.metqtl.txt')
met_file = gsub('SNP', snp, met_file)
met = read.table(met_file, sep = '\t', header = T, stringsAsFactors = F)

ld_file = gsub('SNP', snp, 'results/ld/SNP.ld')
ld = read.table(ld_file, sep = '', header = T, stringsAsFactors = F)

# Merge the two files together based on rsIDs, then split them out so they both
# have matching chr/pos:
tmp = left_join(gwas, met, by = c('SNP' = 'rsids')) %>% drop_na()

gwas_plot = tmp %>% select(CHR, BP, SNP, P)
met_plot = tmp %>% select(CHR, BP, SNP, pval)

colnames(met_plot)[4] = 'P'

out_name = gsub('SNP', snp, 'results/locuszoom/COM/SEX/SNP_lz.jpg')
out_name = gsub('COM', compound, out_name)
out_name = gsub('SEX', cohort, out_name)

# Generate LZ
plot_dat = list(GWAS = gwas_plot, METSIM = met_plot)
locus.zoom(data = plot_dat,
		   snp = snp,
		   offset_bp = 500000,
		   genes.data = genes,
		   ld.file = ld,
		   file.name = out_name,
		   nplots = TRUE)

