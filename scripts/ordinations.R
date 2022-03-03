library(tidyverse)
library(phyloseq)
library(vegan)
library(funfuns)
library(cowplot)




OTU <- phyloseq::import_mothur(mothur_shared_file = './data/FS12b.shared') %>% t()

TAX <- phyloseq::import_mothur(mothur_constaxonomy_file = './data/FS12b_OTU.taxonomy')
colnames(TAX) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family' , 'Genus' )

MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  mutate(treatment=factor(treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))) %>% 
  column_to_rownames(var='sample_ID') %>% 
  sample_data()

FS12b <- phyloseq(MET, TAX, OTU)

###

FS12b <- prune_taxa(taxa_sums(FS12b) > 5, FS12b)

min(sample_sums(FS12b))
hist(sample_sums(FS12b), breaks=100)


### Ordinations ###
# not including ordinations 
# including ordinations because a reviewer asked for them
# fecal only ordinations #
FS12b_feces <- FS12b %>%
  prune_samples(samples = FS12b@sam_data$tissue =='F')

FS12b_feces_meta <- FS12b_feces@sam_data %>% data.frame()
FS12b_feces_OTU <- rrarefy(FS12b_feces@otu_table, min(rowSums(FS12b_feces@otu_table))) %>%
  data.frame()

FS12b_feces_meta <- 
  FS12b_feces_meta %>%
  mutate(treatment=fct_recode(treatment,CON='Control', RPS='RPS', FAM='Acid', RCS='RCS'), 
         set=paste(day, treatment, sep='_'))


FS12b_feces_nmds <- NMDS_ellipse(metadata = FS12b_feces_meta,
                                 OTU_table = FS12b_feces_OTU,
                                 grouping_set = 'set',distance_method = 'bray')

x <- envfit(ord=FS12b_feces_nmds[[3]], env=FS12b_feces_nmds[[1]]$log_sal)
plot(FS12b_feces_nmds[[3]])
plot(x)

centroid_dat <-
  FS12b_feces_nmds[[1]] %>%
  group_by(day, treatment) %>% 
  summarise(MDS1=mean(MDS1), MDS2=mean(MDS2)) %>% 
  ungroup() %>% 
  mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% 
  arrange(day)




library(ggrepel)

FS12b_feces_nmds[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=treatment))+
  geom_point(data = subset(FS12b_feces_nmds[[1]], select = -treatment), col = "grey")+
  geom_point(alpha=.5) +
  geom_segment(aes(xend=centroidX, yend=centroidY), alpha=.5)+
  geom_path(data=centroid_dat, aes(group=treatment), size=1.2) + 
  geom_label(data=centroid_dat,size=2.2,
             aes(label=sub('D', '',day), fill=treatment),
             color='white', 
             fontface='bold') + 
  theme_bw() + 
  facet_wrap(~treatment, ncol=1) +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  +
  labs(caption = paste('NMDS stress =', signif(FS12b_feces_nmds[[3]]$stress, 2)))+
  theme(legend.position='bottom')


ggsave('output/spectrum_revisions/figure_S1_time_NMDS.jpeg', 
       height=9, width = 7, units = 'in', bg='white')


df_ell <-
  FS12b_feces_nmds[[2]] %>% 
  mutate(MDS1=NMDS1, MDS2=NMDS2, 
         day=gsub('(.*)_(.*)','\\1', group), 
         treatment=gsub('(.*)_(.*)','\\2', group)) %>% 
  mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% 
  arrange(da)
  

# df_ell$day <- gsub('(.*)_(.*)','\\1', df_ell$group)
# df_ell$treatment <- gsub('(.*)_(.*)','\\2', df_ell$group)

FS12b_feces_nmds[[1]] %>%
  mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% 
  arrange(day) %>% 
  ggplot(aes(x=MDS1, y=MDS2, color=treatment))+
  geom_point(data = subset(FS12b_feces_nmds[[1]], select = -day), col = "grey")+
  geom_point(alpha=.5, size=2) +
  geom_segment(aes(xend=centroidX, yend=centroidY), alpha=.5)+
  geom_path(data=df_ell, aes(group=group), size=1.1) +
  theme_bw() + 
  facet_wrap(~day, ncol=1) +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  +
  labs(caption = paste('NMDS stress =', signif(FS12b_feces_nmds[[3]]$stress, 2)))+
  theme(legend.position='bottom')
  

FS12b_feces_nmds[[3]]$stress



ggsave('output/spectrum_revisions/figure_S2_treatment_NMDS.jpeg', 
       height=9, width = 6, units = 'in', bg='white')
