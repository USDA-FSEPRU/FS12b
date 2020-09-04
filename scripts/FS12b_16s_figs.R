## libraries ##
library(tidyverse)
library(phyloseq)
library(vegan)
library(funfuns)
library(broom)
library(DESeq2)

## Build phyloseq object ##

OTU <- read_tsv('./data/FS12b_OTU.tsv') %>%
  column_to_rownames(var='sample_ID') %>%
  as.matrix() %>% 
  otu_table(taxa_are_rows = FALSE)


TAX <- read_tsv('./data/FS12b_tax.tsv') %>%
  column_to_rownames(var = 'OTU') %>%
  as.matrix() %>% 
  tax_table()


MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  column_to_rownames(var='sample_ID') %>% 
  sample_data()

FS12b <- phyloseq(MET, TAX, OTU)

###



### PASTE ZONE ###


FS12b_feces <- FS12b %>% prune_samples(samples = FS12b@sam_data$tissue =='F')

FS12b_feces_meta <- FS12b_feces@sam_data %>% data.frame()
FS12b_feces_OTU <- rrarefy(FS12b_feces@otu_table, min(rowSums(FS12b_feces@otu_table))) %>%
  data.frame()


FS12b_feces_nmds <- NMDS_ellipse(metadata = FS12b_feces_meta,
                                 OTU_table = FS12b_feces_OTU,
                                 grouping_set = 'set',distance_method = 'bray')


FS12b_metanmds <- NMDS_ellipse(metadata = FS12b@sam_data,
                               OTU_table = data.frame(rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table)))),
                               grouping_set = 'set',distance_method = 'bray')

FS12b_metanmds




nums <- FS12b_metanmds[[1]] %>% group_by(set) %>% summarise(N=n())
FS12b_metanmds[[1]]
FS12b_metanmds[[3]]


x <- envfit(ord=FS12b_feces_nmds[[3]], env=FS12b_feces_nmds[[1]]$AULC)
# envfit(ord=FS12b_metanmds[[3]], env=FS12b_metanmds[[1]]$log_sal)
plot(FS12b_feces_nmds[[3]])
plot(x)


FS12b_feces_nmds[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=treatment))+
  geom_point() +
  geom_text(aes(label=pignum))


# transform_sample_counts()
# rarefy_even_depth(FS12b)@otu_table

# I DONT THNK THIS DIST IS FROM A RRAREFIED OTU TABLE #
# Fixed #
FS12b_jac <- vegdist(rarefy_even_depth(FS12b)@otu_table, method = 'bray')
FS12b_jac

attr(FS12b_jac, which = 'Labels') == FS12b@sam_data$sample_ID
dispers <- betadisper(FS12b_jac, group = FS12b@sam_data$set)
pdispers <- permutest(dispers, pairwise = TRUE)
# pdispers$pairwise$observed

dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
FS12b@sam_data$ID == dispersdf$group


# meta$sample_ID %in% dispersdf$group

colnames(dispersdf)[2] <- 'ID'
FS12b_meta <- FS12b_metanmds[[1]]

FS12b_meta <- merge(FS12b_meta, dispersdf, by='ID')
FS12b_meta$day <- as.numeric(gsub('D', '', FS12b_meta$day))
FS12b_meta$dayfact <- factor(FS12b_meta$day)


# FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('Control', 'RPS', 'Acid','ZnCu', 'RCS', 'Bglu'))
FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))


FS12b_meta$shan <- diversity(rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))))

FS12b_meta$rich <- specnumber(rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))))
FS12b_meta$even <- FS12b_meta$shan/log(FS12b_meta$rich)


#fecal shannon

FS12b_meta %>% 
  filter(tissue == 'F') %>%
  ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~day, nrow = 1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Fecal Shannon Diversity (alpha) over time')  + 
  geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)


#fecal even
FS12b_meta %>%
  filter(tissue == 'F') %>%
  ggplot(aes(x=treatment, y=even, group=set, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~day, nrow = 1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ggtitle('Fecal evenness over time')+
  geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)

#fecal rich
FS12b_meta %>% 
  filter(tissue == 'F') %>% 
  ggplot(aes(x=treatment, y=rich, group=set, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~day,  nrow = 1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Fecal richness (num OTUs) over time')+
  geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)


#fecal dispersion
FS12b_meta %>% 
  filter(tissue == 'F') %>%
  ggplot(aes(x=treatment, y=dispers.distances, group=set, fill = treatment)) + 
  geom_boxplot() +
  facet_wrap(~day,  nrow = 1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + 
  ggtitle('Fecal community dispersion over time')+
  geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)


# FS12b_meta %>%
#   filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=dispers.distances, group=pignum, color=treatment)) + 
#   geom_line()

# FS12b_meta %>% filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=shan, group=pignum, color=treatment)) + geom_line()

# FS12b_meta %>% filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=dispers.distances, group=pignum, color=treatment)) +
#   geom_line() + geom_point()

# FS12b_meta %>% filter(tissue == 'F') %>%
#   ggplot(aes(x=day, y=rich, group=pignum, color=treatment)) +
#   geom_line() + geom_point()


#
#### CHANGE TO ANOVA TESTS FOR ALPHA AND DISPERSION ####
# 
# get_pairs <- function(df){
#   pp <- pairwise.wilcox.test(df$dispers.distances, df$treatment, p.adjust.method = 'none')
#   ppdf <- as.data.frame(pp$p.value)
#   ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
#   names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
#   ps
#   
# }

disper_fecal_tests <- FS12b_meta %>%
  filter(tissue =='F') %>%
  group_by(day) %>% 
  nest() %>%
  mutate(ANOVA = map(data, ~ aov(data=., formula = dispers.distances ~ treatment)), 
         TUK   = map(ANOVA, TukeyHSD), 
         tid_tuk=map(TUK, tidy)) %>%
  select(day, tid_tuk) %>% unnest(cols = c(tid_tuk))# %>% select(day, starts_with('control'))

disper_fecal_tests %>% filter(grepl('Control', contrast)) %>% 
  filter(adj.p.value < 0.05)


# 
# get_pairs <- function(df){
#   pp <- pairwise.wilcox.test(df$shan, df$treatment, p.adjust.method = 'none')
#   ppdf <- as.data.frame(pp$p.value)
#   ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
#   names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
#   ps
# 
# }

# 
# shan_fecal_tests <- FS12b_meta %>% filter(tissue =='F') %>% group_by(day) %>% 
#   nest() %>% mutate(pps = map(data, get_pairs)) %>%
#   select(day, pps) %>% unnest() %>% select(day, starts_with('control'))


shan_fecal_tests <- FS12b_meta %>%
  filter(tissue =='F') %>%
  group_by(day) %>% 
  nest() %>%
  mutate(ANOVA = map(data, ~ aov(data=., formula = shan ~ treatment)), 
         TUK   = map(ANOVA, TukeyHSD), 
         tid_tuk=map(TUK, tidy)) %>%
  select(day, tid_tuk) %>% unnest(cols = c(tid_tuk))

shan_fecal_tests %>% filter(grepl('Control', contrast)) %>% 
  filter(adj.p.value < 0.05)


###



shan_fecal_tests <- FS12b_meta %>%
  filter(tissue =='F') %>%
  group_by(day) %>% 
  nest() %>%
  mutate(ANOVA = map(data, ~ aov(data=., formula = shan ~ treatment)), 
         TUK   = map(ANOVA, TukeyHSD), 
         tid_tuk=map(TUK, tidy)) %>%
  select(day, tid_tuk) %>% unnest(cols = c(tid_tuk))

shan_fecal_tests %>% filter(grepl('Control', contrast)) %>% 
  mutate(p.FDR=p.adjust(adj.p.value)) %>% 
  filter(adj.p.value < 0.05)

# 
# library(lme4)
# library(lmerTest)
# 
# lmer(data = FS12b_meta, formula = )
# 


# 
# # dispersion tissues
# FS12b_meta %>% filter(tissue == 'X') %>% ggplot(aes(x=treatment, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))
# 
# FS12b_meta %>% filter(tissue == 'C') %>% ggplot(aes(x=treatment, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))
# 
# FS12b_meta %>% filter(tissue == 'I') %>% ggplot(aes(x=treatment, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))
# 
# # shannon tissues
# FS12b_meta %>% filter(tissue == 'X') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))
# 
# FS12b_meta %>% filter(tissue == 'C') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))
# 
# FS12b_meta %>% filter(tissue == 'I') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))
# 
# 



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



######### PW ADON HERE ########

#should split fecal and tissue 1st to reduce # of permutations...

####### TEMP SO I CAN PLAY WITH DIST MEASURES #####

set.seed(71)

# 

# all_pwad <- pairwise.adonis(data.frame(FS12a_rare@otu_table), FS12a_rare@sam_data$set, perm = 999, sim.method = 'bray', binary = FALSE)

# 
# pwad_to_cont <- all_pwad[grep('Control', all_pwad$pairs),]
# # same_day <- pwad_to_cont[grep('.*_(.*)_.*_.* vs .*_\\1_.*_.*', pwad_to_cont$pairs),]
# same_day_tissue <- pwad_to_cont[grep('(.*)_(.*)_.* vs \\1_\\2_.*', pwad_to_cont$pairs),]
# same_day_tissue$treatment <- sub('D[0-9]+_[FX]_([A-Za-z_]+) vs .*_.*_.*', '\\1',same_day_tissue$pairs)
# 
# same_day_tissue[same_day_tissue$treatment == 'Control',]$treatment <- sub('.*_.*_.* vs .*_.*_(.*)','\\1', same_day_tissue[same_day_tissue$treatment == 'Control',]$pairs)
# # sub('(D[0-9]+)_([FX])_([A-Za-z_]+) vs .*_.*_.*', '\\2',same_day_tissue$pairs)
# same_day_tissue$tissue <- sub('(D[0-9]+)_([FX])_([A-Za-z_]+) vs .*_.*_.*', '\\2',same_day_tissue$pairs)
# same_day_tissue$day <- sub('(D[0-9]+)_([FX])_([A-Za-z_]+) vs .*_.*_.*', '\\1',same_day_tissue$pairs)
# 
# same_day_tissue$tissue <- ifelse(same_day_tissue$tissue == 'F', 'feces', 'cecal_mucosa')
# 
# # same_day_tissue$p.adjusted <- p.adjust(same_day_tissue$p.value, method = 'fdr')
# 
# same_day_tissue <- same_day_tissue %>% mutate(set=paste(tissue, day, sep = '_'))
# same_day_tissue <- same_day_tissue %>% group_by(set) %>% mutate(p.adjusted = p.adjust(p.value, method = 'fdr'))
# same_day_tissue$p.plot <- ifelse(same_day_tissue$p.adjusted <= 0.05, paste('p=', round(same_day_tissue$p.adjusted, 3), sep = ''), NA)
# 
# same_day_tissue$set <- factor(same_day_tissue$set, levels = c('feces_D0', 'feces_D23', 'feces_D30', 'cecal_mucosa_D30'))
# 
# #### NEED TO COME UP WITH COLOR TO USE
# # 
# 
# 
# 
# 
# same_day_tissue %>%
#   ggplot(aes(x=treatment, y=F.Model, fill=treatment)) + 
#   geom_col(color='black') + facet_wrap(~set) + geom_text(aes(label=p.plot), nudge_y = .2) + 
#   ggtitle('Community differences compared to controls', subtitle = 'larger F = greater difference.  pvalues shown when <0.05') + 
#   scale_fill_brewer(palette = 'Set1')

# FS12b <- FS12b %>% prune_samples(samples = sample_data(FS12b)$tissue != 'Q')

# FS12b_rare <- FS12b %>% prune_samples(samples = sample_data(FS12b)$tissue != 'I')
# changed this to only include feces
FS12b_rare <- FS12b %>% prune_samples(samples = sample_data(FS12b)$tissue == 'F')
FS12b_rare <- rarefy_even_depth(FS12b_rare)


min(sample_sums(FS12b_rare))
PW.ad <- pairwise.adonis(x=data.frame(FS12b_rare@otu_table), factors = FS12b_rare@sam_data$set, sim.method = 'bray', p.adjust.m = 'none', perm = 999)

# PW.ad <- pairwise.adonis(x=rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))), factors = FS12b@sam_data$set, sim.method = 'jaccard', p.adjust.m = 'none', permutations = 9999)


###### prob doesnt matter... #####

# report this with beginning diffs in beta div
# adonis(data.frame(FS12b_rare@otu_table) ~ tissue + day + treatment, data = data.frame(FS12b_rare@sam_data))

adonis(data.frame(FS12b_rare@otu_table) ~ day + treatment, data = data.frame(FS12b_rare@sam_data))



#######


PW.ad$pairs


goods <- PW.ad[grep('(.*)_(.*)_(.*)_(.*) vs (.*)_\\2_\\3_(.*)', PW.ad$pairs),]

times <- PW.ad[grep('(.*)_(.*)_(.*)_(.*) vs (.*)_(.*)_\\3_\\4', PW.ad$pairs),]

# length(goods[,1])

# goods$p.adjusted <- p.adjust(p=goods$p.value,method = 'holm')

D0 <- goods[grep('D0', goods$pairs),]
D0$day <- 0

D2 <- goods[grep('D2_', goods$pairs),]
D2$day <- 2
D7 <- goods[grep('D7', goods$pairs),]
D7$day <- 7
D14 <- goods[grep('D14', goods$pairs),]
D14$day <- 14
D21 <- goods[grep('D21', goods$pairs),]
D21$day <- 21

fin <- rbind(D0, D2, D7, D14, D21)

fin$pairs <- gsub('X12b_', '', fin$pairs)
fin$pairs <- gsub('_F_', ' feces ', fin$pairs)
fin$pairs <- gsub('_C_', ' cec_cont ', fin$pairs)
fin$pairs <- gsub('_X_', ' cec_muc ', fin$pairs)
fin$pairs <- gsub('_I_', ' il_muc ', fin$pairs)
fin$pairs <- gsub('_Q_', ' tet ', fin$pairs)


# write.csv(fin, 'mothur_PERMANOVA_results.csv')


to_conts <- fin[grep('Control', fin$pairs),]

to_conts$tissue <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\1', to_conts$pairs)

to_conts$treatment <- gsub('D[0-9]+ ([A-Za-z_]+) ([A-Za-z]+) vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\2', to_conts$pairs)

to_conts$treatment[to_conts$treatment == 'Control'] <- gsub('D[0-9]+ ([A-Za-z_]+) ([A-Za-z]+) vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\3', to_conts[to_conts$treatment == 'Control',]$pairs)



# to_conts$p.hoch <- p.adjust(to_conts$p.value, method = 'hoch')
# to_conts$p.holm <- p.adjust(to_conts$p.value, method = 'holm')
to_conts$p.fdr <- p.adjust(to_conts$p.value, method = 'fdr')
to_conts$p.fdr <- round(to_conts$p.fdr, digits = 3)
to_conts$p.fdr.lab <- ifelse(to_conts$p.fdr < 0.05, to_conts$p.fdr, NA)

to_conts$treatment <- factor(to_conts$treatment, levels=c('RPS', 'Acid', 'ZnCu','RCS', 'Bglu'))


######## FIGURE 3 ##########

# ADD ALPHA DIV and DISPERSION 

to_conts %>% filter(tissue == 'feces') %>% ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
  geom_line(size=1.52) + geom_point(shape=21) + scale_color_manual(values=c('#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_label(color='black') +
  scale_fill_manual(values=c('#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ggtitle('Community differences compared to control group over time', subtitle = )


to_conts


# tmp_adon <- same_day_tissue %>% filter(day %in% c('D0', 'D23'))
# 
# tmp_adon <- tmp_adon[,colnames(tmp_adon) %in% colnames(to_conts)]
# 
# tmp_adon2 <- to_conts[,colnames(to_conts) %in% colnames(tmp_adon)]
# 
# long_adon <- rbind(tmp_adon, tmp_adon2)
# 
# long_adon$day <- sub('D0',-30,long_adon$day)
# long_adon$day <- sub('D23',-7,long_adon$day)
# long_adon$day <- as.numeric(long_adon$day)
# 
# long_adon$p.fdr <- p.adjust(long_adon$p.value, method = 'fdr')
# long_adon$p.fdr <- round(long_adon$p.fdr, digits = 3)
# long_adon$p.fdr.lab <- ifelse(long_adon$p.fdr < 0.05, long_adon$p.fdr, NA)
# 
# long_adon$treatment <- factor(long_adon$treatment, levels=c('RPS', 'Acid', 'ZnCu','RCS', 'Bglu'))
# 
# 
# long_adon %>% filter(tissue == 'feces' & treatment %in% c('RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu')) %>% ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
#   geom_line(size=1.52) + geom_point(shape=21) + scale_color_manual(values=c('#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_label(color='black') +
#   scale_fill_manual(values=c('#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   ggtitle('Community differences compared to control group over time', subtitle = )
# 
# 


#SHOULD PLOT DIVERSITY AND RICHNESS AT SAME TIMEPOINTS HERE



########### HIGH LOW TO CONT ##########

##SHOULD ORDINATE THIS TOO##

FS12b_HL <- FS12b %>% subset_samples(treatment %in% c('Control', 'RPS') & tissue =='F')
FS12b_HL <- FS12b %>% subset_samples(treatment %in% c('Control', 'RPS'))

# FS12b_HL %>% subset_samples(treatment == 'RPS') %>% sample_data() %>% select(pignum)


FS12b_HL@sam_data$shed <- ifelse(FS12b_HL@sam_data$pignum %in% c(373,321,181,392,97), 'low', 
                                 ifelse(FS12b_HL@sam_data$pignum %in% c(50, 93,355, 244), 'high', 'Control'))

FS12b_HL@sam_data$set <- paste(FS12b_HL@sam_data$day, FS12b_HL@sam_data$tissue, FS12b_HL@sam_data$shed, sep = '_')

FS12b_HL@sam_data


PW.ad <- pairwise.adonis(x=rrarefy(FS12b_HL@otu_table, min(rowSums(FS12b_HL@otu_table))), factors = FS12b_HL@sam_data$set, sim.method = 'bray', p.adjust.m = 'none', perm = 9999)
# PW.ad <- pairwise.adonis(x=rrarefy(FS12b_HL@otu_table, min(rowSums(FS12b_HL@otu_table))), factors = FS12b_HL@sam_data$set, sim.method = 'jaccard', p.adjust.m = 'none', permutations = 9999)
PW.ad$pairs


goods <- PW.ad[grep('(.*)_(.*) vs \\1_.*', PW.ad$pairs),]

# times <- PW.ad[grep('(.*)_(.*)_(.*)_(.*) vs (.*)_(.*)_\\3_\\4', PW.ad$pairs),]

length(goods[,1])

# goods$p.adjusted <- p.adjust(p=goods$p.value,method = 'holm')

goods$pairs


D0 <- goods[grep('D0', goods$pairs),]
D0$day <- 0

D2 <- goods[grep('D2_', goods$pairs),]
D2$day <- 2
D7 <- goods[grep('D7', goods$pairs),]
D7$day <- 7
D14 <- goods[grep('D14', goods$pairs),]
D14$day <- 14
D21 <- goods[grep('D21', goods$pairs),]
D21$day <- 21

# I think fin and goods are the same thing right now.... why did I do this again?
fin <- rbind(D0, D2, D7, D14, D21)

fin$pairs <- gsub('X12b_', '', fin$pairs)
fin$pairs <- gsub('_F_', ' feces ', fin$pairs)
fin$pairs <- gsub('_C_', ' cec_cont ', fin$pairs)
fin$pairs <- gsub('_X_', ' cec_muc ', fin$pairs)
fin$pairs <- gsub('_I_', ' il_muc ', fin$pairs)
fin$pairs <- gsub('_Q_', ' tet ', fin$pairs)


# write.csv(fin, 'mothur_PERMANOVA_results.csv')

#within tissues
fin <- fin[grep('.* (.*) .* vs .* \\1 .*', fin$pairs),]


to_conts <- fin[grep('Control', fin$pairs),]

not_conts <- fin[-grep('Control', fin$pairs),]

to_conts$tissue <- gsub('D[0-9]+ (.*) ([A-Za-z_]+) vs D[0-9]+ .* ([A-Za-z]+)', '\\1', to_conts$pairs)
to_conts$treatment <- gsub('D[0-9]+ .* ([A-Za-z_]+) vs D[0-9]+ .* ([A-Za-z]+)', '\\2', to_conts$pairs)


to_conts$p.fdr <- p.adjust(to_conts$p.value, method = 'fdr')
to_conts$p.fdr <- round(to_conts$p.fdr, digits = 3)
to_conts$p.fdr.lab <- ifelse(to_conts$p.fdr < 0.05, to_conts$p.fdr, NA)

# to_conts$treatment <- factor(to_conts$treatment, levels=c('RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu'))

p1 <- to_conts %>% filter(tissue =='feces') %>%  ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
  geom_line(size=1.52) + geom_point(shape=21) + scale_color_brewer(palette = 'Set1') + 
  geom_label(color='black') +
  scale_fill_brewer(palette = 'Set1') + 
  ggtitle('Community differences compared to control group over time', subtitle = 'RPS only') + labs(fill='Shedding', 
                                                                                                     color='Shedding')

p1

to_conts %>% filter(!(tissue %in% c('feces', 'tet'))) %>% 
  ggplot(aes(x=tissue, y=F.Model, fill=treatment)) +
  geom_col(position = 'dodge', color='black') + geom_text(aes(label=p.fdr.lab), position=position_dodge(width = 1), vjust=1.5) + 
  ggtitle('PERMANOVA F.stat. : Difference compared to controls across tissues',
          subtitle = 'Higher values represent a greater difference compared to control')  + scale_fill_brewer(palette = 'Set1')



###

#HIGH LOW ORDINATE#
### need to dephyloseqize these objects before NMDS works
HIGH_LOW_OTU <- rarefy_even_depth(FS12b_HL)@otu_table %>% data.frame()
HIGH_LOW_META <- FS12b_HL@sam_data %>% data.frame()

HIGH_LOW_NMDS <- NMDS_ellipse(OTU_table=HIGH_LOW_OTU, metadata = HIGH_LOW_META, grouping_set = 'set')

HIGH_LOW_NMDS[[1]]$shed <- factor(HIGH_LOW_NMDS[[1]]$shed, levels = c('high', 'low', 'Control'))

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D0') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 0, RPS high/low & control')


HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D2') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 2, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D7') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 7, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D14') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 14, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'X' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Cecal mucosa')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'C' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Cecal contents')

p <- HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'I' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Ileal mucosa')



ggplot2::ggplot_build(p)
# ########### groups compared to their D0  #######
# I think I'll cut this... 
# T0s <- times[grep('D0', times$pairs),]
# 
# T0s$pairs <- gsub('X12b_', '', T0s$pairs)
# T0s$pairs <- gsub('_F_', ' feces ', T0s$pairs)
# 
# #4DAF4A, #377EB8)
# 
# #377EB8
# 
# #E41A1C
# 
# 
# 
# T0s$tissue <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\1', T0s$pairs)
# T0s$treatment <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\2', T0s$pairs)
# 
# T0s$day <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs (D[0-9]+) [A-Za-z_]+ ([A-Za-z]+)', '\\2', T0s$pairs)
# 
# # what's this for??
# T0s[T0s$day == "D0",]$day <- gsub('(D[0-9]+) ([A-Za-z_]+) [A-Za-z]+ vs (D[0-9]+) [A-Za-z_]+ ([A-Za-z]+)', '\\1', T0s[T0s$day == "D0",]$pairs)
# 
# T0s$day <- factor(gsub('D','',T0s$day), levels = c(2,7,14,21))
# T0s$pairs
# 
# T0s$p.fdr <- round(p.adjust(T0s$p.value, 'fdr'),3)
# T0s$p.fdr.lab <- ifelse(T0s$p.fdr <0.05, T0s$p.fdr, NA)
# T0s$treatment <- factor(T0s$treatment, levels = c('Control', 'RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu'))
# 
# 
# # this seems different.... investigate
# 
# T0s %>% filter(tissue == 'feces') %>% ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
#   geom_line(size=1.52) + geom_point(shape=21) + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_label(color='black') +
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+ 
#   ggtitle("Community differences compared to each group's Day 0 conformation", 
#           subtitle = 'FDR corrected pvalues shown in boxes') + xlab('Day (vs Day 0)')
# 










################################ RPS SPLIT #########################
FS12b@sam_data$pignum
FS12_RPS <- subset_samples(FS12b, treatment == 'RPS')
FS12_RPS@sam_data$day <- factor(FS12_RPS@sam_data$day, levels = c('D0', 'D2','D7', 'D14', 'D21'))

FS12_RPS@sam_data$shed <- ifelse(FS12_RPS@sam_data$pignum %in% c(373,321,181,392,97), 'low', 'high')
FS12_RPS@sam_data$shed <- factor(FS12_RPS@sam_data$shed, levels = c('high', 'low'))


FS12_RPS@sam_data$set <- paste(FS12_RPS@sam_data$set, FS12_RPS@sam_data$shed, sep = '_')
#

# Keeping samples separate by day #

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D0' & tissue == 'F')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')

# when is the resultsNames not returning what I expect now?
# shed_low_vs_high is what i expect but now 'shed1' is what is returned...
# ok now its fine, must have ran some code twice?
resultsNames(FS12.de)
tmpres <- results(FS12.de, name = 'shed_low_vs_high', cooksCutoff = FALSE)
tmpres <- lfcShrink(FS12.de, res=tmpres, coef = 'shed_low_vs_high', type = 'apeglm')
tmpres[tmpres$padj < 0.1,]


D0_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                              phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                              name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)

### D0 Q
# 
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D0' & tissue == 'Q')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# D0_Q_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                 phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
#                                 name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)
# 

##### D2
FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D2')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D2_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                             phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                             name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)

#### D7
FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D7')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D7_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                             phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                             name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)


# D14 #

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D14')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D14_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                               phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                               name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)

##### D21 F

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'F')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21F_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                                phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                                name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)



#### Tissue X


FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'X')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21X_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                               phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                               name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)



##### tissue C



FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'C')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21C_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                               phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                               name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)

##### tissue I

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'I')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21I_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                                phyloseq.object = FS12b.glom, pvalue = .1, alpha = 0.1,
                                name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)





#
D0_highlow[[2]]$set <- 'D0_feces'
D2_highlow[[2]]$set <- 'D2_feces'
D7_highlow[[2]]$set <- 'D7_feces'
D14_highlow[[2]]$set <- 'D14_feces'
D21F_highlow[[2]]$set <- 'D21_feces'
D21C_highlow[[2]]$set <- 'D21_cecal_content'
D21X_highlow[[2]]$set <- 'D21_cecal_mucosa'
D21I_highlow[[2]]$set <- 'D21_ileal_mucosa'
#

class(D0_highlow[[2]])


RPS_split_master <- bind_rows(list(D0_highlow[[2]],
                                   D2_highlow[[2]],
                                   D7_highlow[[2]],
                                   D14_highlow[[2]],
                                   D21F_highlow[[2]],
                                   D21C_highlow[[2]],
                                   D21X_highlow[[2]], 
                                   D21I_highlow[[2]]))



RPS_split_master$imp <- ifelse(RPS_split_master$padj <= 0.05, TRUE, FALSE)

RPS_split_master$set <- factor(RPS_split_master$set, levels = c('D0_feces','D2_feces' ,'D7_feces', 'D14_feces', 'D21_feces', 'D21_cecal_content', 'D21_cecal_mucosa', 'D21_ileal_mucosa'))
RPS_split_master <- RPS_split_master %>% mutate(newp2=paste0('p=', newp))
# devtools::install_github('jtrachsel/ggscinames')
library(ggscinames)
library(grid)

RPS_split_master %>% filter(set %in% c('D0_feces' ,'D7_feces', 'D14_feces', 'D21_feces')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1')

RPS_split_master %>% filter(set %in% c('D21_feces', 'D21_cecal_content', 'D21_cecal_mucosa')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1') + xlab('') + labs(fill='Shedding')

RPS_split_master %>% filter(set %in% c('D21_ileal_mucosa')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1') + xlab('') + labs(fill='Shedding')




library(cowplot)

p <- RPS_split_master %>%
  group_by(OTU, Treatment) %>%
  filter(padj <= 0.05) %>%  tally() %>%
  ggplot(aes(x=OTU, y=n, fill=Treatment)) + geom_col() +
  scale_fill_brewer(palette = 'Pastel1') + ylab('occurences') + ggtitle('Number of times OTUs are significantly enriched (p<0.05)\n in either shedding phenotype') + 
  theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('')


RPS_split_master %>%
  group_by(OTU, Treatment) %>%
  tally() %>%
  ggplot(aes(x=OTU, y=n, fill=Treatment)) + geom_col() +
  scale_fill_brewer(palette = 'Pastel1') + ylab('occurences') + ggtitle('Number of times OTUs are significantly enriched (p<0.05)\n in either shedding phenotype') + 
  theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('')


ggplot2::ggplot_build(p)



c("#E41A1C", "#FBB4AE", "#377EB8", "#B3CDE3")



# # MERGE THIS WITH TAX AND PRINT TABLE
# RPS_split_master %>%
#   group_by(OTU, Treatment) %>%
#   filter(padj <= 0.1) %>%  tally()


RPS_split_master <- RPS_split_master %>% mutate(p2=ifelse(padj <= 0.05, 'p < 0.05',
                                                          ifelse(padj <= 0.1, 'p < 0.1', NA)))
RPS_split_master <- RPS_split_master %>% mutate(group=factor(paste(p2, Treatment), levels=c("p < 0.05 high", 
                                                                                            "p < 0.1 high", 
                                                                                            "p < 0.05 low", 
                                                                                            "p < 0.1 low")))


RPS_split_master %>% group_by(OTU, group) %>% tally() %>% 
  ggplot(aes(x=OTU, y=n, fill=group)) + geom_col(color='black') +
  scale_fill_manual(values = c("#E41A1C", "#FBB4AE", "#377EB8", "#B3CDE3")) +
  ylab('occurences') +
  ggtitle('Number of times OTUs are enriched \n in either RPS shedding phenotype') + 
  theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('') + coord_flip()



int_OTUs <- RPS_split_master %>% group_by(OTU, group) %>% tally() %>% filter(n>1) %>% select(OTU) %>% unlist(use.names = FALSE)

# write_csv(RPS_split_ints, 'RPS_split_int_OTUs.csv')
RPS_split_ints <- RPS_split_master %>% filter(OTU %in% int_OTUs) %>% 
  select(OTU, Treatment, Genus) %>% unique()

tax <- as.data.frame(FS12b.glom@tax_table)
tax$OTU <- rownames(tax)




#############################################
# regular Differential abundance 


# should make a function for this....
# it would take timepoint, tissue, and return the sig diff OTUs in a dataframe
# need to add tissue and timepoint to dataframe before return


# unique(FS12b@sam_data$pignum)
# 
# FS12b@sam_data$treatment
# DESeq

# DESeq_difabund <- function(phyloseq, day, tissue, scientific = TRUE, shrink_type='normal', 
#                            alpha=0.1, cooks_cut=FALSE, pAdjustMethod='BH'){
#   
#   # FS12b.glom <- tax_glom(FS12b, taxrank = 'Genus')
#   FS12b.glom <- prune_samples(x = phyloseq, samples = phyloseq@sam_data$day == day & phyloseq@sam_data$tissue == tissue)
#   FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
#   FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~treatment)
#   FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
#   
#   finres <- list()
#   resind <- 1
#   for (i in 2:length(resultsNames(FS12.de))){
#     print(resultsNames(FS12.de)[i])
#     treat <- sub('treatment_(.*)_vs_Control','\\1',resultsNames(FS12.de)[i])
#     comp <- sub('treatment_', '', resultsNames(FS12.de)[i])
#     
#     # i dont think these two strategies for results calc are compatible....
#     res <- results(object = FS12.de, name = resultsNames(FS12.de)[i], alpha=alpha, cooksCutoff = cooks_cut, pAdjustMethod = pAdjustMethod)
#     res <- lfcShrink(FS12.de, coef = resultsNames(FS12.de)[i], type = shrink_type)
#     sigtab = res[which(res$padj < alpha), ]
#     
#     if (nrow(sigtab) != 0){
#       # browser()
#       sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FS12b.glom)[rownames(sigtab), ], "matrix"))
#       sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = scientific)
#       sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
#       sigtab$OTU <- rownames(sigtab)
#       sigtab$tissue <- tissue
#       sigtab$day <- day
#       sigtab$comp <- comp
#       finres[[resind]] <- sigtab
#       
#       resind <- resind + 1
#     }
#     
#     
#     
#   }
#   
#   finres <- bind_rows(finres)
#   return(finres)
#   
# }

# something is afoot......
# apparently I removed the Q tissues.....



tocont <- list(DESeq_difabund(phyloseq = FS12b, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               # DESeq_difabund(phyloseq = FS12b, day = 'D0', tissue = 'Q', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D2', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D7', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D14', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'C', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'X', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'I', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'))


tocont <- bind_rows(tocont)

tocontf <- tocont[abs(tocont$log2FoldChange) > .75,]

tocontf %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + geom_hline(yintercept = 20, color='red', size=3)


#### ON TO SOMETHIGN HERE ####
# some variation of this figure for  dif ab.
# maybe do one panel for fecal difabund
# one panel for tissue difabund
tocontf %>% ggplot(aes(x=Genus, y=log2FoldChange, color=Treatment, shape=tissue)) + 
  geom_point() + coord_flip() + geom_hline(yintercept = 0, color='black', size=1) +
  facet_wrap(~day, scales = 'free')


###
tocontf %>% filter(tissue == 'F') %>% ggplot(aes(x=Genus, y=log2FoldChange, color=Treatment)) + 
  geom_point() + coord_flip() + geom_hline(yintercept = 0, color='black', size=1) +
  facet_wrap(~day, scales = 'free' ,ncol = 5)






biguns <- tocontf %>% group_by(OTU) %>% summarise(tot=sum(log2FoldChange)) %>% filter(tot >20) %>% select(OTU) %>% unlist()

tocont %>% filter(OTU %in% biguns) %>% select(OTU,Genus) %>% unique()

tocontf %>% group_by(OTU, Treatment) %>% tally() %>% filter(n>2) %>% as.data.frame()


tocontf %>% group_by(comp) %>% tally() %>% as.data.frame() %>%
  ggplot(aes(x=comp, y=n)) + geom_col() + ggtitle('number of differentially abundant OTUs over the entire experiment')


#### Ok that wasnt so bad.
# Now, which OTUs changed at Salmonella infection?
# I think I need to add sum_sal info at beginning....

# c(c('D0', 'D2'))

#### log_sal as continuous covariate #### 
formula(paste('~', 'log_sal'))
# FS12b@sam_data$day %in% c(day) & FS12b@sam_data$tissue == 'F'
##### SHOULD REALLY LOOK INTO INTERACTION WITH TREATMENT HERE!!!!!!!! 
### OR SUBSET EACH TREATMENT
blarg <- function(phyloseq_obj, day, tissue, covariate, shrink_type='apeglm'){
  form <- formula(paste('~', covariate))
  # print(form)
  FS12b.glom <- phyloseq_obj %>% prune_samples(samples = phyloseq_obj@sam_data$day %in% c(day) & phyloseq_obj@sam_data$tissue == tissue)
  FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
  
  # FS12b.glom@sam_data$log_sal
  
  FS12b.de <- phyloseq_to_deseq2(FS12b.glom, form)
  FS12b.de <- DESeq(FS12b.de, test = 'Wald', fitType = 'parametric')
  
  # these are not both possible.  Right now only lfcshrink is doing anytihng
  res <- results(FS12b.de, cooksCutoff = FALSE, name = covariate)
  res <- lfcShrink(FS12b.de, coef = covariate, type = shrink_type)
  
  # resultsNames(FS12b.de)
  
  res <- res[!is.na(res$padj),]
  res <- res[res$padj < 0.1,]
  sigtab <- res[abs(res$log2FoldChange) > .1 ,]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = TRUE)
  # sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
  sigtab$OTU <- rownames(sigtab)
  sigtab[[covariate]] <- ifelse(sigtab$log2FoldChange >0 , 'increased', 'decreased')
  # sigtab$salm <- ifelse(sigtab$log2FoldChange >0 , 'increased', 'decreased')
  sigtab <- sigtab[order(sigtab$log2FoldChange),]
  sigtab$OTU <- factor(sigtab$OTU, levels = sigtab$OTU)
  sigtab$day <- day
  sigtab$tissue <- tissue
  
  
  p <- sigtab %>% ggplot(aes_string(x='OTU', y='log2FoldChange', fill=covariate)) +
    geom_col(color='black') + coord_flip() + geom_text(aes(label=Genus, y=0))
  
  return(list(p, sigtab))
  
  
}

blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'log_sal')
# nnnn[['test']] <- 'TEST'
# blarg()
# across all treatments

# blarg(phyloseq_obj = FS12b, day = c('D2', 'D7', 'D14', 'D21'), tissue = 'F', covariate = 'log_sal')

# THESE ARE THE ONES AT DAY 2 THAT HAVE A LINEAR RELATIONSHIP WITH SALMONELLA
# LOG2FOLD CHANGE HERE MEANS whut?


# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 2, FS12b.glom)

# FS12b %>% subset_samples(treatment =='RPS') %>% prune_taxa(taxa_sums() > 2)



global_sal_OTUs <- list(blarg(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'log_sal')[[2]], 
                        blarg(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'log_sal')[[2]],
                        blarg(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'log_sal')[[2]],
                        blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'log_sal')[[2]],
                        blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'log_sal')[[2]])

global_sal_OTUs <- bind_rows(global_sal_OTUs)


## this is the filtered OTUs that differ by treatment
tocontf

# These are the OTUs with some kind of linear relationship with log_sal in their tissue/timepoint

increase <- global_sal_OTUs %>% filter(log2FoldChange > 0) # OTUs associated with more salmonella
decrease <- global_sal_OTUs %>% filter(log2FoldChange < 0) # OTUs associated with less salmonella

# THESE l2fc values represent enrichment relative to control
treat_n_sal_increase <- tocontf[tocontf$OTU %in% increase$OTU,] # these are the OTUs that are associated with a treatment and also increased sal
treat_n_sal_decrease <- tocontf[tocontf$OTU %in% decrease$OTU,] # these are the OTUs that are associated with a treatment and also decreased sal


# 
# treat_n_sal_increase <- treat_n_sal_increase %>% mutate(OTU_day_tis=paste(OTU,day, tissue, sep = '_'))
# treat_n_sal_decrease <- treat_n_sal_decrease %>% mutate(OTU_day_tis=paste(OTU,day, tissue, sep = '_'))
# ## THIS GETS TRICKY BECAUSE THE SPECIFIC DAY/TISSUE COMBINATION THESE OTUs ARE ASSOCIATED WITH SALMONELLA DONT NECESSARILY LINE UP WITH 
# WHEN THEY ARE ENRICHED IN TREATMENTS....

# THESE l2fc values represent relationship with salmonella


# THESE  TWO BLOCKS ONLY SHOW WHEN GLOBAL ASSOCIATION WITH SAL AND TREATMENT ENRICHMENT MATCH UP AT SAME TIMEPOINT/TISSUE
global_increase_treat_match <- increase[increase$OTU %in% tocontf$OTU,] %>%
  select(OTU, log2FoldChange, day, tissue) %>%
  mutate(OTU_day_tis=paste(OTU, day, tissue, sep='_'), 
         sal_rel=log2FoldChange) %>%
  select(-log2FoldChange) %>% 
  right_join(treat_n_sal_increase) %>% na.omit()# these are the OTUs that are associated with an increase in sal and also a treatment



global_decrease_treat_match <- decrease[decrease$OTU %in% tocontf$OTU,] %>%
  select(OTU, log2FoldChange, day, tissue) %>%
  mutate(OTU_day_tis=paste(OTU, day, tissue, sep='_'), 
         sal_rel=log2FoldChange) %>%
  select(-log2FoldChange) %>% 
  right_join(treat_n_sal_decrease) %>% na.omit()

big_glob_decrease_treatmatch <- global_decrease_treat_match %>% group_by(OTU) %>% tally() %>% filter(n>1) %>% select(OTU) %>% unlist()

rbind(global_increase_treat_match, global_decrease_treat_match) %>% 
  ggplot(aes(x=OTU, y=sal_rel, fill=Treatment, color=day)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(y=0, x=OTU, label=Genus), color='black') + 
  ylim(-7, 7) + ggtitle('OTUs with a linear relationship to log_sal \n and enriched in any one treatment')


################ 

decrease[decrease$OTU %in% tocontf$OTU,] %>% select(OTU, log2FoldChange, day, tissue)


global_sal_OTUs[!(global_sal_OTUs$OTU %in% tocontf$OTU),] # these are the OTUs that are not associated with a treatment but have an association with log_sal at some timepoint/tissue
big_globs <- global_sal_OTUs[!(global_sal_OTUs$OTU %in% tocontf$OTU),] %>% group_by(OTU) %>% tally() %>% filter(n>1) %>% select(OTU) %>% unlist()

global_sal_OTUs$day <- factor(global_sal_OTUs$day, levels = c('D2', 'D7', 'D14', 'D21'))

# THIS ONE IS OTUS THAT HAVE SIG LIN RELATIONSHIP WITH LOG_SAL at more than 1 time
# no enrich in any treatment relative to control
# THIS ONE!
global_sal_OTUs %>% filter(OTU %in% big_globs) %>%
  ggplot(aes(x=OTU, y=log2FoldChange, fill=day)) +
  geom_hline(yintercept = 0) +
  geom_col(position = position_dodge2(preserve='single')) +
  geom_text_sciname(aes(x = OTU, y=0, sci=Genus), alpha=.5, size=5) +
  coord_flip() + ylim(-3,3) + ggtitle('OTUs with significant linear relationships with log_sal at more than 1 timepoint\n but not associated with any treatment', 
                                      'Log2FoldChange is magnitude of association with salmonella')

# global_sal_OTUs %>% ggplot(aes(x=OTU, y=log2FoldChange, color))

increase %>% group_by(OTU) %>% tally() %>% filter(n>1)
decrease %>% group_by(OTU) %>% tally() %>% filter(n>1)



# blarg(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'log_sal')[[2]]
# # LOOK FOR OTUs associated with treatment that seem to help salmonella status and OTUS associated with treatment that do not 
# blarg(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'log_sal')
# 
# FS12b@sam_data$log_sal
# 
# # blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'X', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'log_sal')
# # blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'I', covariate = 'log_sal')
# 
# # 

# these make less sense.  THey look at linear relationships between OTUs and AULC
# blarg(phyloseq_obj = FS12b, day = 'D0', tissue = 'F', covariate = 'AULC')
# blarg(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'AULC')
# blarg(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'AULC')
# blarg(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'AULC')
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'AULC')
# 
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'X', covariate = 'AULC')
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'AULC')
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'I', covariate = 'AULC')
# 
# 

#### global VFA
# blarg(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'butyrate')
# FS12b_vfa_prune <- prune_samples(x = FS12b , samples = !(FS12b@sam_data$pignum %in% c(6, 265,453, 458, 461, 469, 472)))
# blarg(phyloseq_obj = FS12b_vfa_prune, day = 'D21', tissue = 'C', covariate = 'butyrate')
# blarg(phyloseq_obj = FS12b_vfa_prune, day = 'D21', tissue = 'C', covariate = 'caproate')
# blarg(phyloseq_obj = FS12b_vfa_prune, day = 'D21', tissue = 'C', covariate = 'valerate')
# 
# blarg(phyloseq_obj = FS12b_vfa_prune, day = 'D21', tissue = 'X', covariate = 'butyrate')
# blarg(phyloseq_obj = FS12b_vfa_prune, day = 'D21', tissue = 'X', covariate = 'caproate')
# blarg(phyloseq_obj = FS12b_vfa_prune, day = 'D21', tissue = 'X', covariate = 'valerate')
# 
# 


###
####MOVED FROM ABOVE ####
#### RPS ONLY BLARG ####
### BLARG BY TREATMENT ####
FS12_RPS <- subset_samples(FS12b, treatment == 'RPS')
FS12_control <- subset_samples(FS12b, treatment == 'Control')
FS12_Acid <- subset_samples(FS12b, treatment == 'Acid')
FS12_RCS <- subset_samples(FS12b, treatment == 'RCS')

# FS12_RPS@sam_data$pignum

### CONTROL
# blarg(phyloseq_obj = FS12_control, day = 'D7',tissue = 'F', covariate = 'log_sal')
control_blarg <- bind_rows(list(blarg(phyloseq_obj = FS12_control, day = 'D2',tissue = 'F', covariate = 'log_sal')[[2]],
                                blarg(phyloseq_obj = FS12_control, day = 'D14',tissue = 'F', covariate = 'log_sal')[[2]],
                                blarg(phyloseq_obj = FS12_control, day = 'D21',tissue = 'F', covariate = 'log_sal')[[2]],
                                blarg(phyloseq_obj = FS12_control, day = 'D21',tissue = 'X', covariate = 'log_sal')[[2]]))
# blarg(phyloseq_obj = FS12_control, day = 'D21',tissue = 'C', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_control, day = 'D21',tissue = 'I', covariate = 'log_sal')
control_blarg$treatment <- 'Control'
##### RPS

RPS_blarg <- bind_rows(list(blarg(phyloseq_obj = FS12_RPS, day = 'D2',tissue = 'F', covariate = 'log_sal')[[2]],
                            blarg(phyloseq_obj = FS12_RPS, day = 'D7',tissue = 'F', covariate = 'log_sal')[[2]],
                            blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'X', covariate = 'log_sal')[[2]]))
# blarg(phyloseq_obj = FS12_RPS, day = 'D14',tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'C', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'I', covariate = 'log_sal')
RPS_blarg$treatment <- 'RPS'

RPS_blarg
# tocontf[tocontf[grep('RPS', tocontf$comp),]
tocontf_RPS <- tocontf[grep('RPS', tocontf$comp),]


##### ACID

# blarg(phyloseq_obj = FS12_Acid, day = 'D2',tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_Acid, day = 'D7',tissue = 'F', covariate = 'log_sal')
acid_blarg <- bind_rows(list(blarg(phyloseq_obj = FS12_Acid, day = 'D14',tissue = 'F', covariate = 'log_sal')[[2]],
                             blarg(phyloseq_obj = FS12_Acid, day = 'D21',tissue = 'X', covariate = 'log_sal')[[2]]))
# blarg(phyloseq_obj = FS12_Acid, day = 'D21',tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_Acid, day = 'D21',tissue = 'C', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_Acid, day = 'D21',tissue = 'I', covariate = 'log_sal')
acid_blarg$treatment <- 'Acid'
#### RCS

RCS_blarg <- bind_rows(list(blarg(phyloseq_obj = FS12_RCS, day = 'D2',tissue = 'F', covariate = 'log_sal')[[2]],
                            blarg(phyloseq_obj = FS12_RCS, day = 'D7',tissue = 'F', covariate = 'log_sal')[[2]],
                            blarg(phyloseq_obj = FS12_RCS, day = 'D21',tissue = 'C', covariate = 'log_sal')[[2]]))
# blarg(phyloseq_obj = FS12_RCS, day = 'D21',tissue = 'X', covariate = 'log_sal')))
# blarg(phyloseq_obj = FS12_RCS, day = 'D14',tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_RCS, day = 'D21',tissue = 'F', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_RCS, day = 'D21',tissue = 'I', covariate = 'log_sal')
RCS_blarg$treatment <- 'RCS'

master_blarg <- rbind(control_blarg, RPS_blarg, acid_blarg, RCS_blarg)

treat_blarg_bigs <- master_blarg %>% group_by(OTU) %>% tally() %>% filter(n>1) %>% select(OTU) %>% unlist()

master_blarg %>% filter(OTU %in% treat_blarg_bigs & abs(log2FoldChange) > .25 & tissue == 'F') %>%
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treatment)) +
  geom_col(color='black') + geom_text(aes(label=Genus, y=0), color='black') + coord_flip() +
  ggtitle('Fecal OTUs with linear relationships to Salmonella within treatment groups')

master_blarg[master_blarg$OTU %in% tocontf$OTU,] %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=treatment)) +
  geom_col(color='black') + geom_text(aes(label=Genus, y=0), color='black') + coord_flip() +
  ggtitle('OTUs with linear relationships to Salmonella within treatment groups \n and significant enrichment in one group relative to control', 
          'LFC values represent relationship with salmonella')


# do this one except only include otus with negative lin rel to sal
# maybe scale size to mimick abs lin rel to sal?
tocontf[tocontf$OTU %in% master_blarg$OTU,] %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_point(color='black', shape=21) + geom_text(aes(label=Genus, y=0), color='black') + coord_flip() +# ylim(-20, 60) +
  ggtitle('OTUs significantly enriched treatment groups \nthat also have a significant linear relationship with salmonella', 
          'LFC indicates enrichment relative to control')


p1 <- blarg(phyloseq_obj = FS12_RPS, day = 'D2',tissue = 'F', covariate = 'log_sal')[[1]]
p2 <- blarg(phyloseq_obj = FS12_RPS, day = 'D7',tissue = 'F', covariate = 'log_sal')[[1]]
p3 <- blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'X', covariate = 'log_sal')[[1]]
# p1 <- blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'X', covariate = 'log_sal')[[1]]
# p1 <- blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'X', covariate = 'log_sal')[[1]]

p1 + ggtitle('OTUs with linear relationship to salmonella \nRPS group, D2 Feces')
p2 + ggtitle('OTUs with linear relationship to salmonella \nRPS group, D7 Feces')
# THIS ONE IS V INTERESTING
p3 + ggtitle('OTUs with linear relationship to salmonella \nRPS group, D21 Cecal tissue')

##### I think this is now a repeat?
### THIS SECTION CALCULATES ALL THE OTUs IN THE RPS GROUP THAT HAVE A LINEAR ASSOCIATION WITH salmonella
# NEEDS blarg function defined below...

# in the case of the log_sal covariate these are matched 16S and salmonella culturing samples
# that is the log_sal is measured from the exact same tissue that the 16S data comes from
# in the case of AULC, the 16S samples are related back to the one AULC fecal shedding value calculated for each pig

D2_f_RPS_log_sal <- blarg(phyloseq_obj = FS12_RPS, day = 'D2',tissue = 'F', covariate = 'log_sal')
D7_f_RPS_log_sal <- blarg(phyloseq_obj = FS12_RPS, day = 'D7',tissue = 'F', covariate = 'log_sal')
#blarg(phyloseq_obj = FS12_RPS, day = 'D14',tissue = 'F', covariate = 'log_sal')
#blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'F', covariate = 'log_sal')

# blarg(phyloseq_obj = FS12_RPS, day = 'D0',tissue = 'F', covariate = 'AULC')
# blarg(phyloseq_obj = FS12_RPS, day = 'D2',tissue = 'F', covariate = 'AULC')
# blarg(phyloseq_obj = FS12_RPS, day = 'D7',tissue = 'F', covariate = 'AULC')
# D14_f_RPS_AULC <-  blarg(phyloseq_obj = FS12_RPS, day = 'D14',tissue = 'F', covariate = 'AULC')
# D21_f_RPS_AULC <- blarg(phyloseq_obj = FS12_RPS, day = 'D21',tissue = 'F', covariate = 'AULC')

D21_x_RPS_log_sal <- blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='X', covariate = 'log_sal') # interesting....
# blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='C', covariate = 'log_sal')
# blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='I', covariate = 'log_sal')



# blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='X', covariate = 'AULC')
# blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='C', covariate = 'AULC')
# blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='I', covariate = 'AULC')


blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='C', covariate = 'butyrate')
blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='C', covariate = 'valerate')
blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='C', covariate = 'caproate')


meta$butyrate

blarg(phyloseq_obj = FS12_RPS, day = 'D0', tissue='F', covariate = 'butyrate')
blarg(phyloseq_obj = FS12_RPS, day = 'D0', tissue='F', covariate = 'caproate')
blarg(phyloseq_obj = FS12_RPS, day = 'D0', tissue='F', covariate = 'valerate')

blarg(phyloseq_obj = FS12_RPS, day = 'D21', tissue='X', covariate = 'caproate')


#########





FS12b.glom  = transform_sample_counts(FS12b, function(x) x / sum(x) )
FS12b.glom = filter_taxa(FS12b.glom, function(x) mean(x) > 1e-5, TRUE)




# PSMELT AND BOXPLOTS HERE!!!!!!!!!
# prune_taxa()

######### WARNIGN!!!!! CAREFUL HERE !!!!!!
# D14 doesnt work here because all RCS pigs have exactly the same shedding level at D14
# FS12b.glom <- FS12b %>% prune_samples(samples = FS12b@sam_data$day != 'D0' & FS12b@sam_data$tissue =='F')
# 
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 5, FS12b.glom)
# 
# FS12b.glom@sam_data$log_sal
# 
# # FS12b.glom@sam_data$log_sal
# 
# FS12b.de <- phyloseq_to_deseq2(FS12b.glom, ~log_sal)
# FS12b.de <- DESeq(FS12b.de, test = 'Wald', fitType = 'parametric')
# 
# resultsNames(FS12b.de)
# 
# # res <- results(FS12b.de, cooksCutoff = FALSE, name = 'log_sal')
# res <- lfcShrink(FS12b.de, coef = 'log_sal', type = 'apeglm')
# 
# # resultsNames(FS12b.de)
# 
# res <- res[!is.na(res$padj),]
# res <- res[res$padj < 0.05,]
# sigtab <- res[abs(res$log2FoldChange) > .1 ,]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FS12b.glom)[rownames(sigtab), ], "matrix"))
# sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = TRUE)
# # sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
# sigtab$OTU <- rownames(sigtab)
# sigtab$salm <- ifelse(sigtab$log2FoldChange >0 , 'increased', 'decreased')
# sigtab <- sigtab[order(sigtab$log2FoldChange),]
# sigtab$OTU <- factor(sigtab$OTU, levels = sigtab$OTU)
# 
# 
# sigtab %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=salm)) +
#   geom_col(color='black') + coord_flip() + geom_text(aes(label=Genus, y=0))
# 
### END WRAP ###



# sigtab$tissue <- tissue
# sigtab$day <- day
# sigtab$comp <- comp
# finres[[resind]] <- sigtab

# merge(FS12b@sam_data, sum_sal, by='pignum')




# below here might be on to something... LRT stuff

#ALL FECES
# D0 vs D2 within treatments
# D0 vs D7 within treatments
# D0 vs D14 within treatments
# D0 vs D21 within treatments
unique(FS12b@sam_data$pignum)

FS12b@sam_data$day <- factor(FS12b@sam_data$day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))
FS12b@sam_data$pignum <- factor(FS12b@sam_data$pignum)


FS12b.glom <- prune_samples(x = FS12b, samples = FS12b@sam_data$treatment == 'Control' & FS12b@sam_data$tissue == 'F')
FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~pignum + day)
FS12.de <- DESeq(FS12.de, test = 'LRT', reduced = ~ pignum)

resultsNames(FS12.de)

test2 <- results(object = FS12.de, name = 'day_D2_vs_D0')
test7 <- results(object = FS12.de, name = 'day_D7_vs_D0')

sigtab2 <- test2[which(test2$padj < 0.1),]
sigtab7 <- test7[which(test7$padj < 0.1),]
sigtab2$log2FoldChange
sigtab7$log2FoldChange

all(rownames(sigtab2) == rownames(sigtab7))


sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(FS12b.glom)[rownames(sigtab2), ], "matrix"))
sigtab2$newp <- format(round(sigtab2$padj, digits = 3), scientific = TRUE)
sigtab2$Treatment <- ifelse(sigtab2$log2FoldChange >=0, 'Salmonella', 'Control')
sigtab2$OTU <- rownames(sigtab2)
sigtab2$tissue <- 'feces'
sigtab2$day <- 2
# sigtab$comp <- comp


sigtab2


### I THINK IM ON TO SOMETHING HERE.

### IDENTIFY IMPORTANT OTUS THAT SEEM TO CHANGE WITH SAL AND THEN PLOT TIME COURSE INFO
### CAN DO BY TREATMENT BUT ALSO CAN DO HIGH LOW SHEDDER SPLIT

######## SUM SAL DF ########





# ################### pig trips ##############
# 
# 
# min(rowSums(FS12b@otu_table))
# 
# 
# test <- data.frame(FS12b@otu_table)
# rownames(test)
# rowSums(test)
# 
# FS12.otu.rare <- rrarefy(test, min(rowSums(test)))
# 
# 
# 
# bray.dist <- vegdist(FS12.otu.rare, method = 'jaccard', binary = FALSE)
# 
# FS12b_meta$sample_ID == rownames(FS12.otu.rare)
# 
# #FS12b_meta <- data.frame(FS12b@sam_data)
# 
# 
# 
# #FS12b_meta$
# 
# FS12b_meta$shan2 <- diversity(FS12b@otu_table, base = 2)
# # FS12b_meta$shan <- diversity(FS12b@otu_table)
# # FS12b_meta$invsimp <- diversity(FS12b@otu_table, index = 'invsimpson')
# 
# #FS12b_meta$day <- factor(FS12b_meta$day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))
# #FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('control', 'RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu'))
# 
# # FS12b_meta %>% filter(tissue == 'F') %>% 
# #   ggplot(aes(x=day, y=shan2, group=set, fill=treatment)) +
# #   geom_boxplot() + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + geom_text(aes(label=pignum))
# 
# # FS12b_meta %>% filter(tissue == 'F') %>% 
# #   ggplot(aes(x=day, y=invsimp, group=set, fill=treatment)) +
# #   geom_boxplot() + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + geom_text(aes(label=pignum))
# # 
# 
# 
# 
# # FS12b_meta %>% filter(day == 'D21') %>% ggplot(aes(x=tissue, y=shan2, group=set, fill=treatment)) + geom_boxplot() + geom_text(aes(label=pignum))
# 
# 
# # FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=day, y=shan, group=set, fill=treatment)) + geom_boxplot()
# # FS12b_meta %>% filter(day == 'D21') %>% ggplot(aes(x=tissue, y=shan, group=set, fill=treatment)) + geom_boxplot()
# 
# 
# 
# 
# #ggplot(FS12b_meta, aes(x=treatment, y=shan, group=set)) + geom_boxplot()
# 
# # min(rowSums(shareds_test))
# # hist(rowSums(shareds_test))
# # sort(rowSums(shareds_test))
# 
# dist.data <- as.data.frame(as.matrix(bray.dist))
# dist.data$from <- rownames(dist.data)
# 
# dist.gather <- gather(data = dist.data, key = 'to', value = 'distance', -from)
# 
# #
# 
# 
# dist.gather$fromPig <- gsub('([X12]+[ab]?)([NP]+)([0-9]+)([dWXDi]+)([0-9]+)([A-Z]?)', '\\3', dist.gather$from)
# 
# #
# 
# dist.gather$FT <- paste(dist.gather$from, dist.gather$to, sep = ' ')
# 
# 
# 
# #dist.gather$TF <- paste(dist.gather$to, dist.gather$from, sep = ' ')
# 
# ######
# # all pig pairwise #
# 
# total_ground_covered <- dist.gather[grep('X12bP([0-9]+)D[0-9]+F X12bP\\1D[0-9]+F', dist.gather$FT),] %>% group_by(fromPig) %>% summarise(allpw=sum(distance),
#                                                                                                                                          num=n())
# 
# rooms <- read.csv('./data/Rooms.csv')
# total_ground_covered$treatment <- ifelse(total_ground_covered$fromPig %in% rooms$X6, 'control',
#                                          ifelse(total_ground_covered$fromPig %in% rooms$X7, 'RPS', 
#                                                 ifelse(total_ground_covered$fromPig %in% rooms$X8, 'Acid', 
#                                                        ifelse(total_ground_covered$fromPig %in% rooms$X9, 'Zn+Cu',
#                                                               ifelse(total_ground_covered$fromPig %in% rooms$X10, 'RCS',
#                                                                      ifelse(total_ground_covered$fromPig %in% rooms$X11, 'Bglu', 'asdfsa'))))))
# 
# 
# 
# 
# sum_sal$fromPig <- sum_sal$pignum
# total_ground_covered$fromPig
# 
# total_ground_covered <- total_ground_covered %>% filter(num == 25)
# 
# boxplot(total_ground_covered$allpw~total_ground_covered$treatment)
# 
# sum_sal
# total_ground_covered <- merge(total_ground_covered, sum_sal, by = 'fromPig')
# 
# ############### NEED TO READ IN SUM_SAL ################
# ########################################################
# 
# total_ground_covered$treatment.y == total_ground_covered$treatment.x
# total_ground_covered <- total_ground_covered %>% mutate(treatment=treatment.x) %>% select(-treatment.x, -treatment.y)
# 
# #cor.test(total_ground_covered$allpw, total_ground_covered$AULC)
# 
# 
# total_ground_covered %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, allpw, method = 'pearson')$p.value,
#                                                            AULCvTRIP_T=cor.test(AULC, allpw, method = 'pearson')$statistic)
# 
# 
# 
# total_ground_covered %>% 
#   ggplot(aes(x=allpw, y=AULC, fill=treatment, color=treatment)) +
#   geom_point(size=2, shape=21) + geom_smooth(method = 'lm', se=FALSE) +
#   ggtitle('Correlation between cumulative community membership change and cumulative shedding', 
#           subtitle = 'correlation stats: RPS pval = 0.02, control pval = 0.31, Bglu pval = 0.42') +
#   xlab('Cumulative Bray-Curtis distance (presence/abscence)')
# 
# 
# 
# 
# 
# #
# ######
# D0_2 <- dist.gather[grep('X12bP([0-9]+)D0F X12bP\\1D2F', dist.gather$FT),]
# #colnames(D0_2)[1] <- 'sample_ID'
# colnames(D0_2)[3] <- 'D0_2'
# D0_2 <- D0_2[,c(3,4)]
# 
# D2_7 <- dist.gather[grep('X12bP([0-9]+)D2F X12bP\\1D7F', dist.gather$FT),]
# #colnames(D2_7)[1] <- 'sample_ID'
# colnames(D2_7)[3] <- 'D2_7'
# D2_7 <- D2_7[,c(3,4)]
# 
# D7_14 <- dist.gather[grep('X12bP([0-9]+)D7F X12bP\\1D14F', dist.gather$FT),]
# #colnames(D7_14)[1] <- 'sample_ID'
# colnames(D7_14)[3] <- 'D7_14'
# D7_14 <- D7_14[,c(3,4)]
# 
# D14_21 <- dist.gather[grep('X12bP([0-9]+)D14F X12bP\\1D21F', dist.gather$FT),]
# #colnames(D14_21)[1] <- 'sample_ID'
# colnames(D14_21)[3] <- 'D14_21'
# D14_21 <- D14_21[,c(3,4)]
# 
# D0_21 <- dist.gather[grep('X12bP([0-9]+)D0F X12bP\\1D21F', dist.gather$FT),]
# #colnames(D14_21)[1] <- 'sample_ID'
# colnames(D0_21)[3] <- 'D0_21'
# D0_21 <- D0_21[,c(3,4)]
# 
# 
# #full_join(D0_2, D2_7)
# pig_trips <- merge(D0_2, D2_7, all = TRUE, by = 'fromPig')
# pig_trips <- merge(pig_trips, D7_14, all = TRUE, by = 'fromPig')
# pig_trips <- merge(pig_trips, D14_21, all = TRUE, by = 'fromPig')
# pig_trips <- merge(pig_trips, D0_21, all = TRUE, by = 'fromPig')
# 
# #rowSums(pig_trips)
# #pig_trips <- na.omit(pig_trips)
# colnames(pig_trips[,c(2:5)])
# pig_trips$trip <- rowSums(pig_trips[,c(2:5)])
# hist(pig_trips$trip, breaks = 10)
# 
# 
# 
# rooms <- read.csv('../FS12/Rooms.csv')
# 
# # add treatment data.  This probably isn't the best way to do this...
# 
# #colnames(sum_sal)[1] <- 'fromPig'
# 
# library(funfuns)
# 
# #NMDS_ellipse(metadata = meta_test, OTU_table = shareds_test, grouping_set = 'pig_pen')
# ###############################
# 
# # pig_trips$treatment <- ifelse(pig_trips$fromPig %in% rooms$X6, 'control',
# #                               ifelse(pig_trips$fromPig %in% rooms$X7, 'RPS', 
# #                                      ifelse(pig_trips$fromPig %in% rooms$X8, 'Acid', 
# #                                             ifelse(pig_trips$fromPig %in% rooms$X9, 'Zn+Cu',
# #                                                    ifelse(pig_trips$fromPig %in% rooms$X10, 'RCS',
# #                                                           ifelse(pig_trips$fromPig %in% rooms$X11, 'Bglu', 'asdfsa'))))))
# # 
# # 
# 
# pig_trips <- merge(pig_trips, sum_sal, by = 'fromPig')
# 
# boxplot(pig_trips$trip~pig_trips$treatment)
# boxplot(pig_trips$D0_21~pig_trips$treatment)
# 
# pairwise.wilcox.test(x=pig_trips$trip, g=pig_trips$treatment, p.adjust.method = 'none')
# 
# #colnames(sum_sal)[1] <- 'fromPig'
# pig_trips %>% filter(treatment == "RPS") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# 
# pig_trips %>% filter(treatment == "control") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# 
# pig_trips %>% 
#   ggplot(aes(x=treatment, y=trip, fill=treatment)) + 
#   geom_boxplot() +
#   geom_jitter(shape=21, color='black', stroke=1.2, size=2, width = .2) +
#   scale_fill_manual(values = c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   ggtitle('Cumulative change in each individual pigs community stucture over 21 days') + ylab("Cumulative Jaccard distance")
# 
# 
# 
# 
# pig_trips_cor <- pig_trips %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, trip)$p.value,
#                                                 AULCvTRIP_T=cor.test(AULC, trip)$statistic,
#                                                 AULCv02_P=cor.test(AULC, D0_2)$p.value,
#                                                 AULCv02_T=cor.test(AULC, D0_2)$statistic,
#                                                 AULCv27_P=cor.test(AULC, D2_7)$p.value,
#                                                 AULCv27_T=cor.test(AULC, D2_7)$statistic,
#                                                 AULCv714_P=cor.test(AULC, D7_14)$p.value,
#                                                 AULCv714_T=cor.test(AULC, D7_14)$statistic,
#                                                 AULCv1421_P=cor.test(AULC, D14_21)$p.value,
#                                                 AULCv1421_T=cor.test(AULC, D14_21)$statistic)
# 
# 
# 
# # pig_trips %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, D0_21)$p.value,
# #                                                 AULCvTRIP_T=cor.test(AULC, D0_21)$statistic)
# 
# 
# pig_trips %>% filter(treatment == "control") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# pig_trips %>% filter(treatment == "RPS") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# pig_trips %>% filter(treatment == "Bglu") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# pig_trips %>% filter(treatment == "Zn+Cu") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# 
# pig_trips %>% 
#   ggplot(aes(x=trip, y=AULC, fill=treatment, color=treatment)) +
#   geom_point(size=3, shape=21, color='black') + geom_smooth(method = 'lm', se=FALSE, size=2) +
#   ggtitle('Correlation between cumulative community change and cumulative shedding', 
#           subtitle = 'correlation stats: RPS pval = 0.037, control pval = 0.09, Bglu pval = 0.08') +
#   xlab('Cumulative Jaccard distance') + scale_color_manual(values = c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   scale_fill_manual(values = c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
# 
# 
# c('#3399FF', 'orange', 'red', 'grey', 'purple')
# c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')
# 
# #testse <- cor.test(pig_trips$trip, pig_trips$AULC)
# #testse$p.value
# #testse$statistic
# 
# 
# ggplot(pig_trips, aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# ggplot(pig_trips, aes(x=D0_2, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# ggplot(pig_trips, aes(x=D2_7, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# ggplot(pig_trips, aes(x=D7_14, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# ggplot(pig_trips, aes(x=D14_21, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# 
# ggplot(pig_trips_test, aes(x=mean_trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm', fill = NA) + geom_text(aes(label=pignum))
# 
# #ggplot(pig_trips, aes(x=sum, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
# 
# pig_trips %>% group_by(treatment) %>% summarise(num=n())
# 
# apply(X = pig_trips, MARGIN = 2, FUN = mean, na.rm=TRUE)
# 
# mean(pig_trips$D0_2, na.rm=TRUE)
# mean(pig_trips$D2_7, na.rm=TRUE)
# mean(pig_trips$D7_14, na.rm=TRUE)
# mean(pig_trips$D14_21, na.rm=TRUE)
# 
# median(pig_trips$D0_2, na.rm=TRUE)
# median(pig_trips$D2_7, na.rm=TRUE)
# median(pig_trips$D7_14, na.rm=TRUE)
# median(pig_trips$D14_21, na.rm=TRUE)
# 
# 
# # looking for missing samples
# 
# 
# sum(shared_table[grep('P50D0', rownames(shared_table)),])
# sum(shared_table[grep('P181D7F', rownames(shared_table)),])
# 
# 
# 
# ggplot(pig_trips, aes(x=treatment, y=trip, fill=treatment)) +
#   geom_boxplot() + ylab('Cumulative bray-curtis dissimilarity (each pig)') + geom_jitter(size=2.5,width = 0.2, shape=21)+
#   ggtitle('Cumulative change in community structure through Salmonella infection')
# 

########### cor stuff  ###########

# D0 correlations
fec_VFAs <- res.all

fec_VFAs_0 <- fec_VFAs %>% filter(time == 0) %>% mutate(day=time) %>% select(day, everything(),-time)

ttttt <- FS12b_meta %>% group_by(day) %>% nest()

FS12b_meta %>% group_by(day) %>% nest()

colnames(ttttt$data[[1]])





col_nams_map <- function(df){
  colnames(df) <- paste(day)
}
map()


get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$dispers.distances, df$treatment, p.adjust.method = 'none')
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
  
}


shan_fecal_tests <- FS12b_meta %>% filter(tissue =='F') %>% group_by(day) %>% 
  nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(day, pps) %>% unnest() %>% select(day, starts_with('control'))


########## MISSING DATA??  #########



tttt <- FS12b_meta %>%filter(tissue =='F') %>%  group_by(pignum, day) %>% tally() %>% spread(key = day, value = n)

tttt <- FS12b_meta %>% select(pignum, treatment) %>% unique() %>% left_join(tttt, by = 'pignum')

pig_trips %>% ggplot(aes(x=D0_2, y = D2_7)) + geom_point(aes(color = treatment),size=3) + geom_point()
pig_trips %>% ggplot(aes(x=D0_2, y = D7_14)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D0_2, y = D14_21)) + geom_point(aes(color = treatment),size=3)

pig_trips %>% ggplot(aes(x=D2_7, y = D0_2)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D2_7, y = D7_14)) + geom_point(aes(color = treatment),size=3) + geom_smooth(method = 'lm')
pig_trips %>% ggplot(aes(x=D2_7, y = D14_21)) + geom_point(aes(color = treatment),size=3)

pig_trips %>% ggplot(aes(x=D7_14, y = D0_2)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D7_14, y = D2_7)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D7_14, y = D14_21)) + geom_point(aes(color = treatment),size=3)

pig_trips %>% ggplot(aes(x=D14_21, y = D0_2)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D14_21, y = D2_7)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D14_21, y = D7_14)) + geom_point(aes(color = treatment),size=3)


pig_trips$missing <- ifelse(pig_trips$pignum %in% c(50,181,211,240,253,469), TRUE, FALSE)
pig_trips_test <- pig_trips[,c(1:5, 7)]
PTgath <- pig_trips_test %>% gather(key = interval, value = distance, -fromPig)

avtrp <- PTgath %>% group_by(fromPig) %>% summarise(mean_trip = mean(distance, na.rm = TRUE))

pig_trips_test <- merge(pig_trips, avtrp, by = 'fromPig')

pig_trips_test %>% ggplot(aes(x=trip, y=mean_trip)) + geom_point()

#########

phyloseq::transform_sample_counts()


phyloseq::transform_sample_counts()

FS12_RPS <- subset_taxa(FS12_RPS, taxa_sums(FS12_RPS) > 1)

# plot_bar(FS12_RPS, x='shed')


FS12_RPS_sam <- as_data_frame(FS12_RPS@sam_data)

wht <- FS12_RPS_sam %>% group_by(pignum, tissue) %>% tally()
# missing 50 and 181 fecals

FS12_RPS_otu <- as.data.frame(FS12_RPS@otu_table)
FS12_RPS_otu <- FS12_RPS_otu/rowSums(FS12_RPS_otu) # transforms to relative abundance
FS12_RPS_tax <- as.data.frame(FS12_RPS@tax_table)
FS12_RPS_tax$OTU <- rownames(FS12_RPS_tax)
#
colSums(FS12_RPS_otu)
colsums97 <- colSums(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu)),])/nrow(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu)),])
colsums_others <- colSums(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu), invert=TRUE),])/nrow(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu), invert=TRUE),])
lowerin97 <- (colsums97 - colsums_others) < 0
higherin97 <- (colsums97 - colsums_others) > 0
#

FS12_RPS_otu$sample_ID <- rownames(FS12_RPS_otu)



FS12_RPS_all <- merge(FS12_RPS_sam, FS12_RPS_otu, by='sample_ID')

FS12_RPS_all[1:10, 1:10]

FS12_gath <- FS12_RPS_all %>% gather(key=OTU, value=relabund, -(sample_ID:shed))
FS12_RPS_tax


FS12_gath %>% ggplot(aes(x=pignum, y=relabund)) + geom_col()





