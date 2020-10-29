library(tidyverse)
library(phyloseq)


descripts <- read_tsv('./picrust2_rare_out/pathways_out/path_abun_unstrat_descrip.tsv') %>%
  select(pathway, description)



countData <-
  read_tsv('./picrust2_rare_out/pathways_out/path_abun_unstrat.tsv') %>% 
  column_to_rownames(var = 'pathway') %>%
  floor()



otu_tab <- otu_table(countData, taxa_are_rows = TRUE)
metadat <- sample_data(MET)

metadat$treatment <- factor(metadat$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))


picrust_phylo <- phyloseq(otu_tab, metadat)


D0_fec <- 
  picrust_phylo %>%
  prune_samples(samples = picrust_phylo@sam_data$tissue =='C' &
                          picrust_phylo@sam_data$day =='D21') %>% 
  phyloseq_to_deseq2(design = ~ treatment)



D0_deseq <- DESeq(D0_fec)

resultsNames(D0_deseq)

res <- 
  results(D0_deseq, name = "treatment_RPS_vs_Control") %>%
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
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