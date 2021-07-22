library(tidyverse)
library(phyloseq)
library(vegan)
library(funfuns)
library(cowplot)




OTU <- phyloseq::import_mothur(mothur_shared_file = './data/FS12b.shared') %>% t()

TAX <- phyloseq::import_mothur(mothur_constaxonomy_file = './data/FS12b_OTU.taxonomy')
colnames(TAX) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family' , 'Genus' )

# 
# TAX <- 
#   read_tsv('./raw_data/NEW_MOTHUR_OUT/FS12b_OTU.taxonomy') %>%
#   column_to_rownames(var = 'OTU') %>%
#   as.matrix() %>% 
#   tax_table()

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
# fecal only ordinations #
FS12b_feces <- FS12b %>%
  prune_samples(samples = FS12b@sam_data$tissue =='F')

FS12b_feces_meta <- FS12b_feces@sam_data %>% data.frame()
FS12b_feces_OTU <- rrarefy(FS12b_feces@otu_table, min(rowSums(FS12b_feces@otu_table))) %>%
  data.frame()


FS12b_feces_nmds <- NMDS_ellipse(metadata = FS12b_feces_meta,
                                 OTU_table = FS12b_feces_OTU,
                                 grouping_set = 'set',distance_method = 'bray')

# 
# FS12b_metanmds <- NMDS_ellipse(metadata = FS12b@sam_data,
#                                OTU_table = data.frame(rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table)))),
#                                grouping_set = 'set',distance_method = 'bray')

# FS12b_metanmds




# nums <- FS12b_metanmds[[1]] %>% group_by(set) %>% summarise(N=n())
# FS12b_metanmds[[1]]
# FS12b_metanmds[[3]]

x <- envfit(ord=FS12b_feces_nmds[[3]], env=FS12b_feces_nmds[[1]]$log_sal)
# envfit(ord=FS12b_metanmds[[3]], env=FS12b_metanmds[[1]]$log_sal)
plot(FS12b_feces_nmds[[3]])
plot(x)

FS12b_feces_nmds[[1]] %>%
  group_by(day, treatment) %>% 
  summarise(mds1=mean(MDS1), mds2=mean(MDS2)) %>% 
  ungroup() %>% arrange(day) %>% 
  mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% 
  ggplot(aes(x=mds1, y=mds2, color=treatment)) + 
  geom_point() + 
  geom_path(aes(group=treatment)) + 
  geom_text(aes(label=day))

FS12b_feces_nmds[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=treatment))+
  geom_point() +
  geom_text(aes(label=pignum))


df_ell <- FS12b_metanmds[[2]]

df_ell$experiment <- gsub('(.*)_(.*)_(.*)_(.*)','\\1', df_ell$group)
df_ell$day <- gsub('(.*)_(.*)_(.*)_(.*)','\\2', df_ell$group)
df_ell$day <- gsub('D', '', df_ell$day)
df_ell$day <- factor(df_ell$day, levels = c(0, 2, 7, 14, 21))

df_ell$tissue <- gsub('(.*)_(.*)_(.*)_(.*)','\\3', df_ell$group)
df_ell$treatment <- gsub('(.*)_(.*)_(.*)_(.*)','\\4', df_ell$group)



FS12b_meta$day <- factor(FS12b_meta$day, levels = c(0, 2, 7, 14, 21))
FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('Control', 'RPS', 'Acid','ZnCu', 'RCS', 'Bglu'))
df_ell$treatment <- factor(df_ell$treatment, levels = c('Control', 'RPS', 'Acid','ZnCu', 'RCS', 'Bglu'))


greys <- FS12b_meta

### adding feces only coordinates ###

feces_nmds <- FS12b_feces_nmds[[1]]

FS12b_meta <- feces_nmds %>% select(ID, MDS1, MDS2) %>%
  mutate(fMDS1=MDS1, fMDS2=MDS2) %>%
  select(ID, fMDS1, fMDS2) %>%
  right_join(FS12b_meta)


###



FS12b_meta %>% filter(tissue == 'F' & day == 0) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'F' & day == 0), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 0', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 2) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'F' & day == 2), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 2', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 7) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'F' & day == 7), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 7', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 14) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'F' & day == 14), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 14', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 21) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'F' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 21 feces', x=-.5, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'C' & day == 21) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  # geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'C' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 21\nCecal Contents', x=-0, y=.0, size = 7)



FS12b_meta %>% filter(tissue == 'X' & day == 21) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  # geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'X' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 21\nCecal Mucosa', x=-0, y=.6, size = 7)

FS12b_meta %>% filter(tissue == 'I' & day == 21) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  # geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) +
  geom_path(data = filter(df_ell, tissue == 'I' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +
  annotate(geom='text', label='Day 21\nIleal Mucosa', x=-0, y=.6, size = 5)

####
# INSERTED #

### This is interesting... What does this mean?? ####
# need to move this whole section down below merge #
# FS12b_meta %>% filter(tissue == 'F') %>%
#   group_by(day, treatment) %>%
#   summarise(fMDS1=mean(fMDS1),
#             fMDS2=mean(fMDS2)) %>%
#   ggplot(aes(x=day, y=fMDS1, group=treatment, color=treatment)) +
#   geom_line() + geom_point()


##### #####
# FS12b_meta %>% ggplot(aes(pen, fill=treatment))+geom_histogram(binwidth = 1)

# FS12b_meta %>% filter(tissue == 'F') %>%
#   group_by(day, treatment) %>%
#   summarise(fMDS1=mean(fMDS1),
#             fMDS2=mean(fMDS2)) %>%
#   ggplot(aes(x=day, y=fMDS2, group=treatment, color=treatment)) +
#   geom_line() + geom_point()

# FS12b_meta %>% filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=fMDS1, group=pignum)) +
#   geom_line() + geom_point(aes(color=pen))

# FS12b_meta %>% filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=fMDS2, group=pignum, color=treatment)) +
#   geom_line() + geom_point()

# FS12b_meta %>% filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=even, group=pignum, color=treatment)) +
#   geom_line() + geom_point()


