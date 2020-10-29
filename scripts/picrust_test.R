library(tidyverse)
library(phyloseq)

MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  column_to_rownames(var='sample_ID') %>% 
  sample_data()

descripts <- read_tsv('./data/path_abun_unstrat_descrip.tsv') %>%
  select(pathway, description)



countData <-
  read_tsv('./data/path_abun_unstrat.tsv') %>% 
  column_to_rownames(var = 'pathway') %>%
  floor()


countData <- countData[rowSums(countData > 5) > 5,]

otu_tab <- otu_table(countData, taxa_are_rows = TRUE)
metadat <- sample_data(MET)

metadat$treatment <- factor(metadat$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))


picrust_phylo <- phyloseq(otu_tab, metadat)

hist(taxa_sums(picrust_phylo), breaks = 1000)


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

plotCounts(D0_deseq, gene='CENTFERM-PWY', intgroup = 'shed')


res %>% ggplot(aes(x=log2FoldChange, y=pathway)) + geom_col() + 
  geom_text(aes(x=0, label=description)) + xlim(-10, 10)


