

blastcols <- "qseqid sseqid qstart qend sstart send bitscore length pident staxid ssciname stitle"

blastcols <- blastcols %>% strsplit(split = ' ') %>% unlist()

otu_blast <- read_tsv('./data/OTUrep_blast.tsv', col_names = blastcols)


otu_blast_class <- 
  otu_blast %>%
  filter(pident>85) %>% 
  group_by(qseqid) %>%
  slice_max(order_by = bitscore, with_ties=FALSE) %>% 
  ungroup() %>% 
  mutate(title=sub('(\\[?[A-Z][a-z]+\\]? [a-z]+)(.*)?', '\\1', stitle), 
         OTU=qseqid, 
         genus=sub('(.*) (.*)','\\1',title), 
         species=sub('(.*) (.*)','\\2',title)) %>%
  select(OTU, pident, title, genus, species) %>% 
  write_tsv('output/OTU_blast_class.tsv')




otu_blast %>% filter(qseqid == 'Otu00345')



LOOK <- read_tsv('./output/Allstar_OTUs.tsv') %>% left_join(otu_blast_class)

LOOK %>% select(OTU, perc_comm, Phylum:Genus, title, pident) %>% write_tsv('./output/Allstar_OTUs2.tsv')


otu16 <- 
  read_tsv('./data/OTU16.txt',
         comment = '#',
         col_names = c('qacc', 'sacc', 'evalue', 'qstart', 'qend', 'sstart', 'send', 'stitle', 'pident', 'length')) %>% 
  group_by(qacc) %>% 
  slice_max(pident) %>% 
  filter(grepl('Salmonella', stitle))
