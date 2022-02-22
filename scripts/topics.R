library(fastTopics)
library(tidyverse)
library(phyloseq)

###

OTU <- phyloseq::import_mothur(mothur_shared_file = './data/FS12b.shared') %>% t()

TAX <- phyloseq::import_mothur(mothur_constaxonomy_file = './data/FS12b_OTU.taxonomy')
colnames(TAX) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family' , 'Genus' )

MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  mutate(treatment=fct_recode(treatment, CON='Control', RPS='RPS', FAM='Acid', RCS='RCS')) %>% #pull(treatment)
  column_to_rownames(var='sample_ID') %>% 
  sample_data()

FS12b <- phyloseq(MET, TAX, OTU)

###


FS12b_feces <- 
  FS12b %>%
  prune_samples(samples = FS12b@sam_data$tissue =='F' & 
                  FS12b@sam_data$treatment %in% c('CON', 'RPS', 'FAM', 'RCS'))

FS12b_feces <- prune_taxa(taxa_sums(FS12b_feces) > 5, FS12b_feces)
otu_mat <- as(FS12b_feces@otu_table, 'matrix')

library(fastTopics)

k3 <- fit_topic_model(otu_mat, k=3)
k4 <- fit_topic_model(otu_mat, k=4)
k5 <- fit_topic_model(otu_mat, k=5)
k6 <- fit_topic_model(otu_mat, k=6)
k7 <- fit_topic_model(otu_mat, k=7)
k8 <- fit_topic_model(otu_mat, k=8)
k9 <- fit_topic_model(otu_mat, k=9)


fit_topic_model(X = otu_mat)

