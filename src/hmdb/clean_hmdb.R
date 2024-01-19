library(xml2)
library(tidyr)
library(dplyr)
library(purrr)

# Load HMDB XML
dat = read_xml('data/hmdb/hmdb_metabolites.xml')

# Figure out what "names" are available to pull out
xml_chld = xml_children(dat)
xml_var = xml_name(xml_children(xml_chld[1]))

ns = xml_ns_rename(xml_ns(dat), d1 = 'db')

# Variable names of interest:
# - synonyms
# - taxonomy
# - protein associations
# - general references/pubmed ID
# - biological properties

# Make a function to pull out specific data field from the xml, given an xml
# data for a single metabolite
#
# dat = xml data for a single metabolite
# fields = field(s) of interest, listed in order of the xml structure
# ns = xml namespace
# all = if extracting all tag-value pairs in a block, set to TRUE and it will
#       pull out all the values with each column being the tag name
pull_info = function(dat, fields, ns, all = F) {
	search_list = paste('./db:', fields, sep = '')
	input_list = c(list(dat), search_list)
	met_name = xml_find_first(dat, './db:name', ns) %>% xml_text
	res = reduce(input_list, xml_find_all, ns = ns)
	if (all) {
		res = map_dfr(res, ~ data.frame(metabolite = met_name, col_var = xml_name(.x), value = xml_text(.x)))
		if (nrow(res) == 0) return(NULL)
		res = res %>% group_by(col_var) %>% mutate(num = seq(1:n())) %>% ungroup %>% pivot_wider(names_from = col_var, values_from = value) %>% select(!num)
		res = as.data.frame(res)
	} else {
		res = map_dfr(res, ~ data.frame(metabolite = met_name, x = xml_text(.x)))
		colnames(res) = c('metabolite', fields[length(fields)])
	}
	return(res)
}

# First, make a "base" table with all the metabolites. Other information will
# be added to this base table, so we don't miss any metabolites that have
# missing information in the field of interest
base = map_dfr(xml_chld, ~ pull_info(.x, c('name'), ns)) %>% select(metabolite)

# Pull out all the synonyms:
syn = map_dfr(xml_chld, ~ pull_info(.x, c('synonyms', 'synonym'), ns, all = T))
syn = left_join(base, syn) %>% distinct
base_syn = data.frame(metabolite = base$metabolite, synonym = base$metabolite)
syn = rbind(syn, base_syn) %>% distinct

# Save metabolite synonym table:
write.table(syn, 'data/hmdb/hmdb_synonyms.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Pull out the status:
status = map_dfr(xml_chld, ~ pull_info(.x, c('status'), ns))

write.table(status, 'data/hmdb/hmdb_status.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Pull out various chemical "classes":
kdm = map_dfr(xml_chld, ~ pull_info(.x, c('taxonomy', 'kingdom'), ns))
super = map_dfr(xml_chld, ~ pull_info(.x, c('taxonomy', 'super_class'), ns))
cls = map_dfr(xml_chld, ~ pull_info(.x, c('taxonomy', 'class'), ns))
sub = map_dfr(xml_chld, ~ pull_info(.x, c('taxonomy', 'sub_class'), ns))
mol_frm = map_dfr(xml_chld, ~ pull_info(.x, c('taxonomy', 'molecular_framework'), ns))

merged = reduce(list(kdm, super, cls, sub,mol_frm), full_join)

write.table(merged, 'data/hmdb/hmdb_chemtax.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# List of "alternative parents"
alt = map_dfr(xml_chld, ~ pull_info(.x, c('taxonomy', 'alternative_parents', 'alternative_parent'), ns))

write.table(alt, 'data/hmdb/hmdb_alt_parents.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Pull out protein association info:
prot = map_dfr(xml_chld, ~ pull_info(.x, c('protein_associations', 'protein', '*'), ns, all = T))

write.table(prot, 'data/hmdb/hmdb_metprot.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Pull out references
ref = map_dfr(xml_chld, ~ pull_info(.x, c('general_references', 'reference', '*'), ns, all = T))

write.table(prot, 'data/hmdb/hmdb_metprot.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Pull out biological properties/pathways:
cell = map_dfr(xml_chld, ~ pull_info(.x, c('biological_properties', 'cellular_locations', 'cellular'), ns))
write.table(cell, 'data/hmdb/hmdb_cell_loc.txt', sep = '\t', col.names = T, row.names = F, quote = F)

bio = map_dfr(xml_chld, ~ pull_info(.x, c('biological_properties', 'tissue_locations', 'tissue'), ns))
write.table(bio, 'data/hmdb/hmdb_tissue_loc.txt', sep = '\t', col.names = T, row.names = F, quote = F)

pathways = map_dfr(xml_chld, ~ pull_info(.x, c('biological_properties', 'pathways', 'pathway', '*'), ns, all = T))

write.table(pathways, 'data/hmdb/hmdb_biopath.txt', sep = '\t', col.names = T, row.names = F, quote = F)

