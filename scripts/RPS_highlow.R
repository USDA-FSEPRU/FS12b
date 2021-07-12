########## HIGH LOW TO CONT ##########
library(phyloseq)
library(vegan)
library(tidyverse)
library(funfuns)
library(cowplot)


shed_data <- 
  read_tsv('./data/sal_summary.tsv') %>% 
  mutate(shedding=ifelse(pignum %in% c(373,321,181,392,97), 'low', 
         ifelse(pignum %in% c(50, 93,355, 244), 'high', 'control'))) %>% 
  filter(treatment %in% c('control', 'RPS'))

# pbuild <- ggplot_build(p1)

# pbuild$data

my_cols <- c(low='#377EB8', high='#E41A1C', control='#4DAF4A')
#4DAF4A green
#E41A1C red
#377EB8 blue

p1 <- shed_data %>%
  ggplot(aes(x=treatment, y=AULC)) +
  geom_boxplot() +
  geom_jitter(aes(fill=shedding), width=.2, size=3, shape=21) + 
  # scale_fill_brewer(palette = 'Set1') + 
  scale_fill_manual(values = my_cols) + 
  theme_cowplot() + 
  xlab('')



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




##SHOULD ORDINATE THIS TOO##

# FS12b_HL <- FS12b %>% subset_samples(treatment %in% c('Control', 'RPS') & tissue =='F')
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


D0 <- goods[grep('0', goods$pairs),]
D0$day <- 0

D2 <- goods[grep('2_', goods$pairs),]
D2$day <- 2
D7 <- goods[grep('7_', goods$pairs),]
D7$day <- 7
D14 <- goods[grep('14', goods$pairs),]
D14$day <- 14
D21 <- goods[grep('21', goods$pairs),]
D21$day <- 21

# I think fin and goods are the same thing right now.... why did I do this again?
fin <- rbind(D0, D2, D7, D14, D21)

# fin$pairs <- gsub('X12b_', '', fin$pairs)
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
### WRITE TO CONTS ####
# to_conts %>% write_tsv('./fig_dat/High_low_PERMANOVA.tsv')

# figure 5ish
p2 <- to_conts %>% filter(tissue =='feces') %>% 
  ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
  geom_line(size=1.52) +
  geom_point(shape=21, show.legend = F) + 
  scale_color_brewer(palette = 'Set1') + 
  geom_label(color='black', show.legend = F) +
  theme_cowplot()+
  scale_fill_brewer(palette = 'Set1') + 
  labs(color='Shedding') #+
  # xlab('')
  # ggtitle('Community differences compared to control group over time', subtitle = 'RPS only') + labs(fill='Shedding', 
                                                                                                     

p2

p3 <- 
  to_conts %>% filter(!(tissue %in% c('feces', 'tet', 'il_muc'))) %>% 
  ggplot(aes(x=tissue, y=F.Model, fill=treatment)) +
  geom_col(position = 'dodge', color='black') + geom_text(aes(label=p.fdr.lab), position=position_dodge(width = 1), vjust=1.5) + 
  # ggtitle('PERMANOVA F.stat. : Difference compared to controls across tissues',
  #         subtitle = 'Higher values represent a greater difference compared to control')  +
  scale_fill_brewer(palette = 'Set1') + 
  xlab('')+
  theme_cowplot() + 
  labs(fill='Shedding')


p3





to_conts %>% write_tsv('./output/PERMANOVAs_highlow_vs_control.tsv')





### cowplot

fig_S1 <- ggdraw()+
  draw_plot(p1, 0,.45,.45,.55)+
  draw_plot(p3, .45,.45,.55,.45)+
  draw_plot(p2, 0,0,1,.45)+
  draw_plot_label(x=c(0,0, .45), y=c(1,.45,1), label = c('A', 'B','C'))
fig_S1

ggsave('output/figureS1.jpeg', height=5, width = 7, units = 'in')

###

#HIGH LOW ORDINATE#
### need to dephyloseqize these objects before NMDS works
HIGH_LOW_OTU <- rarefy_even_depth(FS12b_HL)@otu_table %>% data.frame()
HIGH_LOW_META <- FS12b_HL@sam_data %>% data.frame()

HIGH_LOW_NMDS <- NMDS_ellipse(OTU_table=HIGH_LOW_OTU, metadata = HIGH_LOW_META, grouping_set = 'set')

HIGH_LOW_NMDS[[1]]$shed <- factor(HIGH_LOW_NMDS[[1]]$shed, levels = c('high', 'low', 'Control'))

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D0') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + 
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 0, RPS high/low & control')


HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D2') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + 
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 2, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D7') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + 
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 7, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D14') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) +
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 14, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) +
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'X' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + 
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Cecal mucosa')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'C' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) +
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Cecal contents')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'I' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + 
  # geom_point(size=3)+
  geom_text(aes(label=pignum))+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Ileal mucosa')






################################ RPS SPLIT #########################


# NEED TO SET SHRINK TYPE CONSISTENTLY
FS12b@sam_data$pignum
FS12_RPS <- subset_samples(FS12b, treatment == 'RPS')
# FS12_RPS@sam_data$day <- factor(FS12_RPS@sam_data$day, levels = c('D0', 'D2','D7', 'D14', 'D21'))

FS12_RPS@sam_data$shed <- ifelse(FS12_RPS@sam_data$pignum %in% c(373,321,181,392,97), 'low', 'high')
FS12_RPS@sam_data$shed <- factor(FS12_RPS@sam_data$shed, levels = c('high', 'low'))


FS12_RPS@sam_data$set <- paste(FS12_RPS@sam_data$set, FS12_RPS@sam_data$shed, sep = '_')
#

# Keeping samples separate by day #

highlow_DESEQ <- 
  function(phyloseq, day, tissue, shrinktype, cookscutoff){
    FS12b.glom <-
      prune_samples(x = phyloseq,
                    samples = phyloseq@sam_data$day == day &
                      phyloseq@sam_data$tissue == tissue)
    
    FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
    FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
    FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
    tmpres <- results(FS12.de, name = 'shed_low_vs_high', cooksCutoff = cookscutoff)
    tmpres <- lfcShrink(FS12.de, res=tmpres, coef = 'shed_low_vs_high', type = shrinktype)
    tmpres$day <- day
    tmpres$tissue <- tissue
    tmpres$OTU <- rownames(tmpres)
    tmpres$Treatment <- ifelse(tmpres$log2FoldChange > 0 , 'low shed', 'high shed')
    
    tmpres = cbind(as(tmpres, "data.frame"), as(tax_table(FS12b.glom)[rownames(tmpres),], "matrix"))
    
    # tmpres <- as(tax_table(phyloseq), Class = 'matrix') %>% as_data_frame() %>%  right_join(tmpres)
    
    return(tmpres)
  }

# D0_highlow <- highlow_DESEQ(FS12_RPS, day = 0, tissue = 'F', shrinktype = 'apeglm', cookscutoff = FALSE)
# 
RPS_split_master <- 
  bind_rows(
    highlow_DESEQ(FS12_RPS, day = 'D0', tissue = 'F', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D2', tissue = 'F', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D7', tissue = 'F', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D14', tissue = 'F', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D21', tissue = 'F', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D21', tissue = 'X', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D21', tissue = 'C', shrinktype = 'apeglm', cookscutoff = TRUE),
    highlow_DESEQ(FS12_RPS, day = 'D21', tissue = 'I', shrinktype = 'apeglm', cookscutoff = TRUE))

# library(IHW)
ihwRes <- IHW::ihw(pvalue ~ baseMean,  data = RPS_split_master, alpha = 0.05)


RPS_split_master$IHW_pval <- IHW::adj_pvalues(ihwRes)

RPS_split_master <- RPS_split_master %>% filter(IHW_pval < 0.1 & abs(log2FoldChange) > 0.5)

#  NEED TO DO WORK HERE #

RPS_split_master$imp <- ifelse(RPS_split_master$IHW_pval <= 0.05, TRUE, FALSE)
RPS_split_master$set <- paste(RPS_split_master$day, RPS_split_master$tissue, sep = '_')
RPS_split_master$set <- factor(RPS_split_master$set, levels = c('D0_F','D2_F' ,'D7_F', 'D14_F', 'D21_F', 'D21_C', 'D21_X', 'D21_I'))
RPS_split_master$newp <- signif(RPS_split_master$IHW_pval, digits = 2)
RPS_split_master <- RPS_split_master %>% mutate(newp2=paste0('p=', newp))
# devtools::install_github('jtrachsel/ggscinames')
library(ggscinames)
library(grid)

# IHW has an alpha thng that messes up ggplot2
# detach("package:IHW", unload=TRUE)

RPS_split_master %>% #filter(set %in% c('D0_feces' ,'D7_feces', 'D14_feces', 'D21_feces')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1')

#### THIS ONE ###
RPS_split_master %>% #filter(set %in% c('D0_feces' ,'D7_feces', 'D14_feces', 'D21_feces')) %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=day, shape=tissue)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21:27))+
  # geom_text_sciname(aes(x=Genus, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) +
  coord_flip() +
  scale_fill_brewer(palette = 'Set1')


# RPS_split_master %>% filter(set %in% c('D21_feces', 'D21_cecal_content', 'D21_cecal_mucosa')) %>%
#   ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
#   geom_bar(stat='identity') + 
#   geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
#   facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1') + xlab('') + labs(fill='Shedding')
# 
# RPS_split_master %>% filter(set %in% c('D21_ileal_mucosa')) %>%
#   ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
#   geom_bar(stat='identity') + 
#   geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
#   facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1') + xlab('') + labs(fill='Shedding')





# p <- RPS_split_master %>%
#   group_by(OTU, Treatment) %>%
#   filter(padj <= 0.05) %>%  tally() %>%
#   ggplot(aes(x=OTU, y=n, fill=Treatment)) + geom_col() +
#   scale_fill_brewer(palette = 'Pastel1') + ylab('occurences') + ggtitle('Number of times OTUs are significantly enriched (p<0.05)\n in either shedding phenotype') + 
#   theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('')
# 
# 
# RPS_split_master %>%
#   group_by(OTU, Treatment) %>%
#   tally() %>%
#   ggplot(aes(x=OTU, y=n, fill=Treatment)) + geom_col() +
#   scale_fill_brewer(palette = 'Pastel1') + ylab('occurences') + ggtitle('Number of times OTUs are significantly enriched (p<0.05)\n in either shedding phenotype') + 
#   theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('')


# ggplot2::ggplot_build(p)



# c("#E41A1C", "#FBB4AE", "#377EB8", "#B3CDE3")



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


# RPS_split_master %>% group_by(OTU, group) %>% tally() %>% 
#   ggplot(aes(x=OTU, y=n, fill=group)) + geom_col(color='black') +
#   scale_fill_manual(values = c("#E41A1C", "#FBB4AE", "#377EB8", "#B3CDE3")) +
#   ylab('occurences') +
#   ggtitle('Number of times OTUs are enriched \n in either RPS shedding phenotype') + 
#   theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('') + coord_flip()


# 
# int_OTUs <- RPS_split_master %>% group_by(OTU, group) %>% tally() %>% filter(n>1) %>% select(OTU) %>% unlist(use.names = FALSE)
# 
# # write_csv(RPS_split_ints, 'RPS_split_int_OTUs.csv')
# RPS_split_ints <- RPS_split_master %>% filter(OTU %in% int_OTUs) %>% 
#   select(OTU, Treatment, Genus) %>% unique()
# 
# tax <- as.data.frame(FS12b.glom@tax_table)
# tax$OTU <- rownames(tax)
# 



#### END HIGHLOW INSERT

### STACKED BARS FOR RPS ###

FS12_RPS@sam_data$pignum[order(FS12_RPS@sam_data$AULC)]

RPS_bars <- FS12_RPS %>% transform_sample_counts(function(x) x/sum(x)) %>% psmelt()


RPS_bars %>%
  mutate(Class_lump= fct_lump_n(Genus, n = 9, w=Abundance)) %>% 
  mutate(pig_fact=fct_reorder(.f = as.factor(pignum), .x = AULC)) %>% 
  ggplot(aes(x=pig_fact, y=Abundance, fill=Class_lump)) +
  geom_col(color=alpha(colour = 'black', alpha = .5)) + facet_wrap(~day+tissue) +
  scale_fill_brewer(palette = 'Set1')



# GPr  = transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
###



