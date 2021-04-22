library(tidyverse)
library(phyloseq)
library(DESeq2)




MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  column_to_rownames(var='sample_ID') %>% 
  sample_data()


# 
# descripts <- read_tsv('./data/no_rare_path_abun_unstrat_descrip.tsv') %>% 
#   select(pathway, description)
# 
# countData <-
#   read_tsv('./data/no_rare_path_abun_unstrat_descrip.tsv') %>% 
#   select(pathway, MET$ID) %>% 
#   column_to_rownames(var = 'pathway') %>%
#   floor()
# 

####


# pathways
descripts <- read_tsv('./data/FS12b_path_descrip.tsv') %>% 
  select(pathway, description)


countData <-
  read_tsv('./data/FS12b_path_descrip.tsv') %>% 
  select(pathway, MET$ID) %>% 
  column_to_rownames(var = 'pathway') %>%
  floor()



# # KO
# descripts <- read_tsv('./data/KO_pred_metagenome_descript.tsv') %>% 
#   transmute(KO=`function`,description=description)
# 
# 
# countData <-
#   read_tsv('./data/KO_pred_metagenome_descript.tsv') %>% 
#   select(`function`, MET$ID) %>% 
#   column_to_rownames(var = 'function') %>%
#   floor()
# 
# 
# # EC
# descripts <- read_tsv('./data/EC_pred_metagenome_unstrat_decrip.tsv') %>% 
#   transmute(KO=`function`,description=description)
# 
# 
# countData <-
#   read_tsv('./data/EC_pred_metagenome_unstrat_decrip.tsv') %>% 
#   select(`function`, MET$ID) %>% 
#   column_to_rownames(var = 'function') %>%
#   floor()
# 



####

hist(rowSums(countData), breaks=100)
hist(rowSums(countData > 5), breaks = 100)
countData <- countData[rowSums(countData > 5) > 5,]
hist(rowSums(countData > 5), breaks = 100)


otu_tab <- otu_table(countData, taxa_are_rows = TRUE)
metadat <- sample_data(MET)
tax_tab <- descripts %>% column_to_rownames('pathway') %>% as.matrix() %>% tax_table()



metadat$treatment <- factor(metadat$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))


picrust_phylo <- phyloseq(otu_tab, metadat, tax_tab)

hist(taxa_sums(picrust_phylo), breaks = 1000)

hist(sample_sums(picrust_phylo), breaks = 1000)

picrust_phylo@sam_data$ID

#
KO_difabund <- 
  function(phyloseq, tissue, day, descripts, pval=0.05, l2fc=.25){
    deseq_obj <- 
      phyloseq %>%
      prune_samples(samples = phyloseq@sam_data$tissue ==tissue &
                      phyloseq@sam_data$day ==day) %>% 
      prune_taxa(taxa = taxa_sums(phyloseq) > 5) %>% 
      phyloseq_to_deseq2(design = ~ treatment) %>%
      DESeq()
    
    res <- 
      lfcShrink(deseq_obj, coef = which(resultsNames(deseq_obj) == 'treatment_RPS_vs_Control')) %>% 
      as.data.frame() %>% 
      filter(padj < pval & abs(log2FoldChange) > l2fc) %>% 
      rownames_to_column(var = 'KO') %>% 
      left_join(descripts) %>% 
      mutate(KO=fct_reorder(KO, log2FoldChange), 
             tissue= tissue, 
             day=day)
    
    which(resultsNames(deseq_obj) == 'treatment_RPS_vs_Control')
    # res <- 
    #   results(deseq_obj, name = "treatment_RPS_vs_Control") %>%
    #   as.data.frame() %>% 
    #   filter(padj < pval & abs(log2FoldChange) > l2fc) %>% 
    #   rownames_to_column(var = 'pathway') %>% 
    #   left_join(descripts) %>% 
    #   mutate(pathway=fct_reorder(pathway, log2FoldChange), 
    #          tissue= tissue, 
    #          day=day)
    
    
    
    p1 <-
      res %>%
      ggplot(aes(x=log2FoldChange, y=KO)) +
      geom_col() +
      geom_text(aes(x=0, label=description)) +
      xlim(-10, 10)
    
    return(list(res, p1))
    
  }


#

picrust_difabund <- 
function(phyloseq, tissue, day, descripts, pval=0.05, l2fc=.25){
  deseq_obj <- 
    phyloseq %>%
    prune_samples(samples = phyloseq@sam_data$tissue ==tissue &
                    phyloseq@sam_data$day ==day) %>% 
    prune_taxa(taxa = taxa_sums(phyloseq) > 5) %>% 
    phyloseq_to_deseq2(design = ~ treatment) %>%
    DESeq()
  
  res <- 
    lfcShrink(deseq_obj, coef = which(resultsNames(deseq_obj) == 'treatment_RPS_vs_Control')) %>% 
    as.data.frame() %>% 
    filter(padj < pval & abs(log2FoldChange) > l2fc) %>% 
    rownames_to_column(var = 'pathway') %>% 
    left_join(descripts) %>% 
    mutate(pathway=fct_reorder(pathway, log2FoldChange), 
           tissue= tissue, 
           day=day)
  
  which(resultsNames(deseq_obj) == 'treatment_RPS_vs_Control')
  # res <- 
  #   results(deseq_obj, name = "treatment_RPS_vs_Control") %>%
  #   as.data.frame() %>% 
  #   filter(padj < pval & abs(log2FoldChange) > l2fc) %>% 
  #   rownames_to_column(var = 'pathway') %>% 
  #   left_join(descripts) %>% 
  #   mutate(pathway=fct_reorder(pathway, log2FoldChange), 
  #          tissue= tissue, 
  #          day=day)
  
  
  
  p1 <-
    res %>%
    ggplot(aes(x=log2FoldChange, y=pathway)) +
    geom_col() +
    geom_text(aes(x=0, label=description)) +
    xlim(-10, 10)
  
  return(list(res, p1))
  
}


RPS_vs_control_picrust <- list(
  picrust_difabund(picrust_phylo, tissue = 'F', day = 'D0', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'F', day = 'D2', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'F', day = 'D7', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'F', day = 'D14', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'F', day = 'D21', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'X', day = 'D21', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'C', day = 'D21', descripts = descripts),
  picrust_difabund(picrust_phylo, tissue = 'I', day = 'D21', descripts = descripts))
# 
# RPS_vs_control_KO <- list(
#   KO_difabund(picrust_phylo, tissue='F', day='D0', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='F', day='D2', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='F', day='D7', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='F', day='D14', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='F', day='D21', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='X', day='D21', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='C', day='D21', descripts = descripts),
#   KO_difabund(picrust_phylo, tissue='I', day='D21', descripts = descripts))
# 
# CONTROL_VS_RPS_KO <- bind_rows(lapply(RPS_vs_control_KO, '[[', 1)) %>% 
#   mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21')), 
#          KO=fct_reorder(KO, log2FoldChange), 
#          description=fct_reorder(description, log2FoldChange))
# 
# 
# CONTROL_VS_RPS_KO %>% filter(day == 'D14') %>% 
#   arrange(desc(log2FoldChange)) %>% 
#   ggplot(aes(x=KO, y=log2FoldChange)) +
#   geom_hline(yintercept = 0)+
#   geom_text(aes(label=description, y=0)) +
#   geom_point(aes(color=day, shape=tissue)) +
#   coord_flip() + 
#   scale_color_brewer(palette = 'Set1') +
#   # facet_wrap(~ptype, scales = 'free_y', ncol = 1) + 
#   theme(panel.grid.major = element_line(color='grey'))
# 
# 



CONTROL_VS_RPS_PICRUST <- bind_rows(lapply(RPS_vs_control_picrust, '[[', 1)) %>% 
  mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21')), 
         pathway=fct_reorder(pathway, log2FoldChange), 
         description=fct_reorder(description, log2FoldChange))%>%
  mutate(ptype=ifelse(grepl('synthesis', description), 'synthesis',
                      ifelse(grepl('degradation', description), 'degradation', 'other')))


# write_tsv(CONTROL_VS_RPS_PICRUST, './output/Control_vs_RPS_picrust.tsv')

CONTROL_VS_RPS_PICRUST %>% filter(ptype == 'other')


CONTROL_VS_RPS_PICRUST %>% #filter(tissue =='F') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  geom_hline(yintercept = 0)+
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~ptype, scales = 'free_y', ncol = 1) + 
  theme(panel.grid.major = element_line(color='grey'))



CONTROL_VS_RPS_PICRUST %>% filter(tissue !='F') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  geom_hline(yintercept = 0)+
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~ptype, scales = 'free_y', ncol = 1) + 
  theme(panel.grid.major = element_line(color='grey'))



library(cowplot)


CONTROL_VS_RPS_PICRUST %>% filter(ptype =='other') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue), size=3) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') +theme_cowplot() + 
  geom_hline(yintercept = 0)
  # facet_wrap(~ptype, scales = 'free_y', ncol = 1)


CONTROL_VS_RPS_PICRUST %>% filter(ptype =='synthesis') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') + theme_cowplot()+
  theme_cowplot() + 
  geom_hline(yintercept = 0)


CONTROL_VS_RPS_PICRUST %>% filter(ptype =='degradation') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue), size=3) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') +theme_cowplot() + 
  geom_hline(yintercept = 0)+theme_cowplot() + 
  geom_hline(yintercept = 0)



# tocont %>% filter(Phylum == 'Proteobacteria')




CONTROL_VS_RPS_PICRUST %>% filter(tissue !='F') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1')



# library(DESeq2)
# 
# pic_path_met <- 
#   PICPATH %>%
#   pivot_longer(names_to = 'ID', cols = starts_with('X')) %>% 
#   pivot_wider(names_from = pathway) %>% left_join(MET)
# 
# 
# pic_path_met
# 




###### highlow cont

# FS12b_HL <- picrust_phylo %>% subset_samples(treatment %in% c('Control', 'RPS') & tissue =='F')
FS12b_HL <- picrust_phylo %>% subset_samples(treatment %in% c('RPS'))

# FS12b_HL %>% subset_samples(treatment == 'RPS') %>% sample_data() %>% select(pignum)


FS12b_HL@sam_data$shed <- ifelse(FS12b_HL@sam_data$pignum %in% c(373,321,181,392,97), 'low', 
                                 ifelse(FS12b_HL@sam_data$pignum %in% c(50, 93,355, 244), 'high', 'Control'))

FS12b_HL@sam_data$set <- paste(FS12b_HL@sam_data$day, FS12b_HL@sam_data$tissue, FS12b_HL@sam_data$shed, sep = '_')


highlow_picrust <- 
  function(phylo_obj, tissue, day, pval = 0.05, descripts){
    deseq_obj <- 
    phylo_obj %>%
    prune_samples(samples = phylo_obj@sam_data$tissue ==tissue &
                    phylo_obj@sam_data$day ==day) %>% 
    prune_taxa(taxa = taxa_sums(phylo_obj) > 5) %>% 
    phyloseq_to_deseq2(design = ~ shed) %>%
    DESeq()
  
  

  # deseq_obj
  
  
  
  # deseq_obj <- DESeq(deseq_obj)
  
  resultsNames(deseq_obj)
  
  
  # lfcShrink(dds = deseq_obj, coef = 2)
  
  res <- 
    lfcShrink(dds = deseq_obj, coef = 2) %>%
    as.data.frame() %>% 
    filter(padj < pval) %>% 
    rownames_to_column(var = 'pathway') %>% 
    left_join(descripts) %>% 
    mutate(pathway=fct_reorder(pathway, log2FoldChange))
  
  print(plotCounts(deseq_obj, gene='PWY-3781', intgroup = 'shed'))
  print(nrow(res))
  if (nrow(res) < 1){
    stop('no sigs')
  }
  
  
  
  p1 <- 
    res %>% ggplot(aes(x=log2FoldChange, y=pathway)) + geom_col() + 
    geom_text(aes(x=0, label=description)) + xlim(-10, 10)
  return(list(p1, res))
}




highlow_KO <- 
  function(phylo_obj, tissue, day, pval = 0.05, descripts){
    deseq_obj <- 
      phylo_obj %>%
      prune_samples(samples = phylo_obj@sam_data$tissue ==tissue &
                      phylo_obj@sam_data$day ==day) %>% 
      prune_taxa(taxa = taxa_sums(phylo_obj) > 5) %>% 
      phyloseq_to_deseq2(design = ~ shed) %>%
      DESeq()
    
    
    
    # deseq_obj
    
    
    
    # deseq_obj <- DESeq(deseq_obj)
    
    resultsNames(deseq_obj)
    
    
    # lfcShrink(dds = deseq_obj, coef = 2)
    
    res <- 
      lfcShrink(dds = deseq_obj, coef = 2) %>%
      as.data.frame() %>% 
      filter(padj < pval & abs(log2FoldChange) > 1) %>% 
      rownames_to_column(var = 'KO') %>% 
      left_join(descripts) %>% 
      mutate(KO=fct_reorder(KO, log2FoldChange))
    
    # print(plotCounts(deseq_obj, gene='PWY-3781', intgroup = 'shed'))
    print(nrow(res))
    if (nrow(res) < 1){
      stop('no sigs')
    }
    
    
    
    p1 <- 
      res %>% ggplot(aes(x=log2FoldChange, y=KO)) + geom_col() + 
      geom_text(aes(x=0, label=description)) + xlim(-10, 10)
    return(list(p1, res))
  }


# highlow_KO(FS12b_HL, tissue = 'F', day = 'D7', descripts = descripts, pval = 0.05)


highlow_picrust(FS12b_HL, tissue = 'F', day = 'D7', descripts = descripts, pval = 0.05)


read_tsv('./data/FS12b_path_descrip.tsv')

HL_phylo <- 
  FS12b_HL %>%
  prune_samples(samples = FS12b_HL@sam_data$tissue =='F' &
                  FS12b_HL@sam_data$day =='D7') %>% 
  prune_taxa(taxa = taxa_sums(FS12b_HL) > 5) %>% 
  phyloseq_to_deseq2(design = ~ shed)





deseq_obj <- DESeq(HL_phylo)

resultsNames(deseq_obj)


lfcShrink(dds = deseq_obj, coef = 2)

res <- 
  lfcShrink(dds = deseq_obj, coef = 2) %>%
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames_to_column(var = 'pathway') %>% 
  left_join(descripts) %>% 
  mutate(pathway=fct_reorder(pathway, log2FoldChange))

nrow(res)

plotCounts(deseq_obj, gene='CENTFERM-PWY', intgroup = 'shed')


res %>% ggplot(aes(x=log2FoldChange, y=pathway)) + geom_col() + 
  geom_text(aes(x=0, label=description)) + xlim(-10, 10)

#################

# linear relationships with log_sal



HL_phylo <- 
  FS12b_HL %>%
  prune_samples(samples = FS12b_HL@sam_data$tissue =='F' &
                  FS12b_HL@sam_data$day =='D7') %>% 
  prune_taxa(taxa = taxa_sums(FS12b_HL) > 5) %>% 
  phyloseq_to_deseq2(design = ~ scale(AULC))


scale(FS12b_HL@sam_data$AULC)


deseq_obj <- DESeq(HL_phylo)

resultsNames(deseq_obj)


lfcShrink(dds = deseq_obj, coef = 2)

res <- 
  lfcShrink(dds = deseq_obj, coef = 2) %>%
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  # filter(log2FoldChange > 0.1) %>% 
  rownames_to_column(var = 'pathway') %>% 
  left_join(descripts) %>% 
  mutate(pathway=fct_reorder(pathway, log2FoldChange))

nrow(res)

plotCounts(deseq_obj, gene='CENTFERM-PWY', intgroup = 'shed')


res %>% ggplot(aes(x=log2FoldChange, y=pathway)) + geom_col() + 
  geom_text(aes(x=0, label=description)) + xlim(-10, 10)





RES <- 
  list(
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D2',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D7',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D14',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01),
# DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01)
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'X', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'valerate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'caproate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'butyrate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'succinate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'propionate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'lactate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'acetate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'phenylacetate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'isobutyrate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'isovalerate', plot_lab = 'description',l2fc_plot = 0.01),
DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'oxalate', plot_lab = 'description',l2fc_plot = 0.01)
)


huh <- RES %>% map(pluck(2)) %>% bind_rows() %>% filter(padj < 0.05)



biosynth <- huh %>% filter(grepl('synthesis', description))
degrade <- huh %>% filter(grepl('degradation', description))
others <- huh %>% filter(!grepl('synthesis', description) & !grepl('degradation', description))


biosynth %>% group_by(OTU, description) %>% tally() %>% arrange(desc(n))
degrade %>% group_by(OTU, description) %>% tally() %>% arrange(desc(n))
others %>% group_by(OTU, description) %>% tally() %>% arrange(desc(n))







######## all now #######

RES <- 
  list(
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D2',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D7',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D14',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'F', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    # DESeq_cov_asso(phyloseq_obj = FS12b_HL, day = 'D21',tissue = 'C', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01)
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'X', covariate = 'log_sal', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'valerate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'caproate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'butyrate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'succinate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'propionate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'lactate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'acetate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'phenylacetate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'isobutyrate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'isovalerate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T),
    DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21',tissue = 'C', covariate = 'oxalate', plot_lab = 'description',l2fc_plot = 0.01, treatment = T)
  )
