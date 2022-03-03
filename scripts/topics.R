# install.packages('topicmodels')
# remotes::install_github("stephenslab/fastTopics")
# devtools::install_github("lasy/alto")


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
  mutate(treatment=fct_recode(treatment, CON='Control', RPS='RPS', FAM='Acid', RCS='RCS'), 
         day_num=as.numeric(sub('D', '', day)),
         day=fct_relevel(day,c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% #pull(treatment)
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


var_params <- 
  list(K1=list(k=1),
       K2=list(k=2),
       K3=list(k=3), 
       K4=list(k=4), 
       K5=list(k=5), 
       K6=list(k=6), 
       K7=list(k=7), 
       K8=list(k=8), 
       K9=list(k=9), 
       K10=list(k=10))

library(alto)
MODS <- alto::run_lda_models(data = otu_mat, lda_varying_params_lists = var_params, dir = 'output/', verbose = T)
ALN <- align_topics(MODS)

plot(ALN)
ggplot(topics(ALN), aes(m, coherence)) +
  geom_point(alpha = 0.5)


META <- MODS$K4$gamma |> as.data.frame() |> rownames_to_column(var='ID') |> left_join(FS12b_feces@sam_data)

# topic 3 seems salmonella related
META |> ggplot(aes(x=V1, y=log_sal, color=treatment)) + geom_point() + facet_wrap(~day)
META |> ggplot(aes(x=V2, y=log_sal, color=treatment)) + geom_point() + facet_wrap(~day)
META |> ggplot(aes(x=V3, y=log_sal, color=treatment)) + geom_point() + facet_wrap(~day)
META |> ggplot(aes(x=V4, y=log_sal, color=treatment)) + geom_point() + facet_wrap(~day)


META |> ggplot(aes(x=day, y=V3, fill=treatment)) + geom_point(shape=21) + geom_boxplot()
