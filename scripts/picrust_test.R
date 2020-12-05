library(tidyverse)
library(phyloseq)
library(DESeq2)




MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  column_to_rownames(var='sample_ID') %>% 
  sample_data()

descripts <- read_tsv('./data/no_rare_path_abun_unstrat_descrip.tsv') %>% 
  select(pathway, description)



countData <-
  read_tsv('./data/no_rare_path_abun_unstrat_descrip.tsv') %>% 
  select(pathway, MET$ID) %>% 
  column_to_rownames(var = 'pathway') %>%
  floor()

hist(rowSums(countData), breaks=100)
hist(rowSums(countData > 5), breaks = 100)
countData <- countData[rowSums(countData > 5) > 5,]
hist(rowSums(countData > 5), breaks = 100)


otu_tab <- otu_table(countData, taxa_are_rows = TRUE)
metadat <- sample_data(MET)

metadat$treatment <- factor(metadat$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))


picrust_phylo <- phyloseq(otu_tab, metadat)

hist(taxa_sums(picrust_phylo), breaks = 1000)

hist(sample_sums(picrust_phylo), breaks = 1000)

picrust_phylo@sam_data$ID

picrust_difabund <- 
function(phyloseq, tissue, day, descripts, pval=0.05, l2fc=.5){
  deseq_obj <- 
    phyloseq %>%
    prune_samples(samples = phyloseq@sam_data$tissue ==tissue &
                    phyloseq@sam_data$day ==day) %>% 
    prune_taxa(taxa = taxa_sums(phyloseq) > 5) %>% 
    phyloseq_to_deseq2(design = ~ treatment) %>%
    DESeq()
  
  
  
  res <- 
    results(deseq_obj, name = "treatment_RPS_vs_Control") %>%
    as.data.frame() %>% 
    filter(padj < pval & abs(log2FoldChange) > l2fc) %>% 
    rownames_to_column(var = 'pathway') %>% 
    left_join(descripts) %>% 
    mutate(pathway=fct_reorder(pathway, log2FoldChange), 
           tissue= tissue, 
           day=day)
  
  
  
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



CONTROL_VS_RPS_PICRUST <- bind_rows(lapply(RPS_vs_control_picrust, '[[', 1)) %>% 
  mutate(day=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21')), 
         pathway=fct_reorder(pathway, log2FoldChange), 
         description=fct_reorder(description, log2FoldChange))%>%
  mutate(ptype=ifelse(grepl('synthesis', description), 'synthesis',
                      ifelse(grepl('degradation', description), 'degradation', 'other')))


write_tsv(CONTROL_VS_RPS_PICRUST, './output/Control_vs_RPS_picrust.tsv')

CONTROL_VS_RPS_PICRUST %>% filter(ptype == 'other')


CONTROL_VS_RPS_PICRUST %>% #filter(tissue =='F') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~ptype, scales = 'free_y', ncol = 1)

CONTROL_VS_RPS_PICRUST %>% filter(ptype =='degradation') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue), size=3) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') #+
  # facet_wrap(~ptype, scales = 'free_y', ncol = 1)


CONTROL_VS_RPS_PICRUST %>% filter(ptype =='synthesis') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') #+


CONTROL_VS_RPS_PICRUST %>% filter(ptype =='degradation') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1') #+



# tocont %>% filter(Phylum == 'Proteobacteria')




CONTROL_VS_RPS_PICRUST %>% filter(tissue !='F') %>% 
  arrange(desc(log2FoldChange)) %>% 
  ggplot(aes(x=description, y=log2FoldChange)) +
  # geom_text(aes(label=description, y=0)) +
  geom_point(aes(color=day, shape=tissue)) + coord_flip() + 
  scale_color_brewer(palette = 'Set1')


D0_fec <- 
  picrust_phylo %>%
  prune_samples(samples = picrust_phylo@sam_data$tissue =='X' &
                          picrust_phylo@sam_data$day =='D21') %>% 
  prune_taxa(taxa = taxa_sums(picrust_phylo) > 5) %>% 
  phyloseq_to_deseq2(design = ~ treatment)



D0_deseq <- DESeq(D0_fec)

resultsNames(D0_deseq)

res <- 
  results(D0_deseq, name = "treatment_RPS_vs_Control") %>%
  as.data.frame() %>% 
  filter(padj < 0.1) %>% 
  rownames_to_column(var = 'pathway') %>% 
  left_join(descripts) %>% 
  mutate(pathway=fct_reorder(pathway, log2FoldChange))



res %>% ggplot(aes(x=log2FoldChange, y=pathway)) + geom_col() + 
  geom_text(aes(x=0, label=description)) + xlim(-10, 10)


library(DESeq2)

pic_path_met <- 
  PICPATH %>%
  pivot_longer(names_to = 'ID', cols = starts_with('X')) %>% 
  pivot_wider(names_from = pathway) %>% left_join(MET)


pic_path_met





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



highlow_picrust(FS12b_HL, tissue = 'F', day = 'D21', descripts = descripts, pval = 0.05)




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
  prune_samples(samples = FS12b_HL@sam_data$tissue =='X' &
                  FS12b_HL@sam_data$day =='D21') %>% 
  prune_taxa(taxa = taxa_sums(FS12b_HL) > 5) %>% 
  phyloseq_to_deseq2(design = ~ log_sal)





deseq_obj <- DESeq(HL_phylo)

resultsNames(deseq_obj)


lfcShrink(dds = deseq_obj, coef = 2)

res <- 
  lfcShrink(dds = deseq_obj, coef = 2) %>%
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  filter(log2FoldChange > 0.1) %>% 
  rownames_to_column(var = 'pathway') %>% 
  left_join(descripts) %>% 
  mutate(pathway=fct_reorder(pathway, log2FoldChange))

nrow(res)

plotCounts(deseq_obj, gene='CENTFERM-PWY', intgroup = 'shed')


res %>% ggplot(aes(x=log2FoldChange, y=pathway)) + geom_col() + 
  geom_text(aes(x=0, label=description)) + xlim(-10, 10)


