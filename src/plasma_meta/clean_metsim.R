library(vroom)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])
dat$cpid = paste(dat$chrom, dat$pos, sep = '_')

chrpos = vroom('data/metsim_data/metsim_var.hg19.vcf', delim = '\t', skip = 98, col_names = c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'))
chrpos$cpid = gsub('chr', '', chrpos$ID)
chrpos$cpid = gsub('_[A-Z]*$', '', chrpos$cpid)
chrpos = chrpos %>% select(CHR, POS, cpid, REF, ALT) %>% mutate(CHR = gsub('chr', '', CHR))

dat = left_join(dat, chrpos, by = c('cpid', 'ref' = 'REF', 'alt' = 'ALT')) %>% select(CHR, POS, ref, alt, pval, beta, sebeta, maf) %>% mutate(cpid = paste(CHR, POS, sep = '_')) %>% filter(!is.na(CHR))
colnames(dat) = c("CHR", "POS", "other_allele", "effect_allele", "P", "effect", "SE", "MAF", "cpid")

# Save output:
out_name = gsub('.tsv.gz', '.clean.txt', args[1])
write.table(dat, out_name, sep = '\t', col.names = T, row.names = F, quote = F)

