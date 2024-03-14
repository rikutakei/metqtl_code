# Genetic colocalisation analysis of Metabolite QTL with gout

## Analysis process

### Colocalisation

1. Take all the European lead gout variants from full, male, and female gout GWAS
1. Based on these variants, identify +/-500kb region (1Mb window) centered around the variant
1. For the full/male/female gout GWAS, pull out all the regions, regardless of where the variants come from (selection of results for the relevant data set will be done later)
1. Repeat above with all of the metabolite GWAS from METSIM and Schlosser data sets
1. Run colocalisation of gout GWAS with the metabolite GWAS, using the subsetted summary stats
1. Merge the colocalisation results

