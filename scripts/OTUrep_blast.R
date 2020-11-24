

blastcols <- "qseqid sseqid qstart qend sstart send bitscore length pident staxid ssciname stitle"

blastcols <- blastcols %>% strsplit(split = ' ') %>% unlist()

otu_blast <- read_tsv('./data/OTUrep_blast.tsv', col_names = blastcols)


otu_blast %>%
  filter(pident>93) %>% 
  group_by(qseqid) %>%
  slice_max(order_by = bitscore, with_ties=FALSE) %>% 
  ungroup() %>% 
  mutate(title=sub('(\\[?[A-Z][a-z]+\\]? [a-z]+)(.*)?', '\\1', stitle), 
         OTU=qseqid, 
         genus=sub('(.*) (.*)','\\1',title), 
         species=sub('(.*) (.*)','\\2',title)) %>%
  select(OTU, pident, title, genus, species) %>% 
  write_tsv('output/OTU_blast_class.tsv')

slice_max()
