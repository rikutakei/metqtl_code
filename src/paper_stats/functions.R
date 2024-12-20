# Function for generating per-locus summary
per_locus = function(dat, metsim = F) {
  # Remove the A/B/C signals from multi-signal loci:
  dat = dat %>% mutate(locus = gsub('_[A-Z]$', '', locus))
  # Merge the SNPs at multi-signal loci:
  dat = dat %>% group_by(compound, phenostring, locus) %>% mutate(SNP = paste(unique(SNP), collapse = ' / ')) %>% ungroup %>% distinct
  # For each locus, count up how many compounds/metabolites colocalised with it
  res = dat %>% distinct(compound, phenostring, SNP, locus) %>% group_by(SNP, locus) %>% summarize(count = n(), compound_id = paste(unique(compound), collapse = '; '), metabolites = paste(unique(phenostring), collapse = '; ')) %>% ungroup %>% arrange(desc(count))
  # It will be cool to see the breakdown of the types of metabolites that
  # colocalised at the locus based on the category information provided in the
  # METSIM data set:
  if (metsim) {
    category_prop = dat %>% arrange(category) %>% mutate(dummy = 1) %>% pivot_wider(names_from = category, values_from = dummy, values_fill = 0) %>% select(compound, phenostring, SNP, locus, Amino_Acid:last_col()) %>% distinct
    category_prop = category_prop %>% group_by(SNP, locus) %>% summarize_at(vars(Amino_Acid:last_col()), sum) %>% ungroup
    res = left_join(res, category_prop)
  }
  return(res)
}

# Simple function to add and assign colours to the categories
add_color = function(dat) {
  categories = c("Amino_Acid", "Carbohydrate", "Cofactors_and_Vitamins", "Energy", "Lipid", "Nucleotide", "Partially_Characterized", "Peptide", "Uncharacterized", "Xenobiotics")
  color = cbind(categories, brewer.pal(10, "Paired")) %>% as.data.frame
  colnames(color) = c('category', 'color')
  res = left_join(dat, color)
  return(res)
}

# Given a plot-ready data set from the per-locus summary table, plot it
plot_perloc = function(dat, facet = F, outname) {
  plot_dat = dat %>% mutate(category = gsub('_', ' ', category))
  plot_dat = plot_dat %>% ggplot(aes(x = locus, y = prop, fill = category))
  if (facet) {
    plot_dat = plot_dat + facet_wrap( ~ category, nrow = 2)
  }
  plot_dat = plot_dat +
    scale_fill_manual(values = unique(dat$color), labels = sort(unique(dat$category)), name = 'Category') +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
    labs(x = '', y = 'Proportion of colocalised metabolite categories at each locus') +
    geom_col()
  ggsave(outname, width = 10, height = 8)
}

# Function to summarise the per-metabolite colocalisation
per_compound = function(dat) {
  # Remove the A/B/C signals from multi-signal loci:
  res = dat %>% mutate(locus = gsub('_[A-Z]$', '', locus))
  # Merge the SNPs at multi-signal loci:
  res = res %>% group_by(compound, phenostring, locus) %>% mutate(SNP = paste(unique(SNP), collapse = '; ')) %>% ungroup %>% distinct
  # Given a unique set of compound-locus pair, count up how many loci
  # colocalised for each metabolite:
  res = res %>% distinct(compound, phenostring, SNP, locus) %>% group_by(compound, phenostring) %>% summarize(count = n(), coloc_variants = paste(unique(SNP), collapse = '; '), loci = paste(unique(locus), collapse = '; ')) %>% ungroup %>% arrange(desc(count))
  return(res)
}

