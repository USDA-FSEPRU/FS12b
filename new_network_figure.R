
library(tidyverse)
library(phyloseq)
library(vegan)
library(funfuns)
library(broom)
library(DESeq2)
library(cowplot)
library(geomnet)

source('./scripts/FS12b_16s_figs.R')

## Build phyloseq object ##
## probably split this into the calculations and then the figures.
### write out the calcs and then read back in in figures script.


# OTU <- 
#   read_tsv('./data/FS12b_OTU.tsv') %>%
#   column_to_rownames(var='Group') %>%
#   as.matrix() %>% 
#   otu_table(taxa_are_rows = FALSE)
OTU <- phyloseq::import_mothur(mothur_shared_file = './raw_data/NEW_MOTHUR_OUT/FS12b.shared') %>% t()

TAX <- phyloseq::import_mothur(mothur_constaxonomy_file = './raw_data/NEW_MOTHUR_OUT/FS12b_OTU.taxonomy')
colnames(TAX) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family' , 'Genus' )

# 
# TAX <- 
#   read_tsv('./raw_data/NEW_MOTHUR_OUT/FS12b_OTU.taxonomy') %>%
#   column_to_rownames(var = 'OTU') %>%
#   as.matrix() %>% 
#   tax_table()

### NEED RPS_sigOTUs

MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  mutate(treatment=factor(treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))) %>% 
  column_to_rownames(var='sample_ID') %>% 
  sample_data()

FS12b <- phyloseq(MET, TAX, OTU)


# split into two functions, one that outputs the network edges/ network object
# one that plots the network?

test_fun <- 
  function(phyloseq_obj, EDGE_QVAL, EDGE_LFC, SEED, highlight_these_otus){
    # browser()
    
    set.seed(SEED)
  FS12b <- phyloseq_obj
  all_melt <- 
    FS12b %>% 
    rarefy_even_depth() %>% 
    psmelt()
  
  taxtab <- as(tax_table(FS12b), 'matrix') %>% data.frame() %>% rownames_to_column(var = 'OTU')
  
  all_agg_abund_info <- 
    all_melt %>%
    group_by(OTU) %>% 
    summarise(mean_abund=mean(Abundance)) %>% 
    mutate(tot_abund=sum(mean_abund), 
           prop_comm=mean_abund/tot_abund, 
           perc_comm=prop_comm * 100) %>% 
    ungroup() %>% 
    left_join(taxtab)
  
  
  # all_agg_abund_info %>% dplyr::count(OTU) %>% filter(n !=1)
  
  
  # padj after all calcs?
  
  sal_OTU_assoc_treat <- 
    list(
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'log_sal', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'X', covariate = 'log_sal', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D0', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'AULC', treatment = T)[[2]]
    ) %>% 
    bind_rows() %>%
    # filter(padj < 0.05 & abs(log2FoldChange) >.5) %>%
    # filter(padj < 0.05 ) %>% 
    mutate(node_name=paste(day, tissue, direction, covariate,sep = '_'))
  
  # sal_OTU_assoc_treat %>% filter(covariate == 'AULC') %>% 
  #   ggplot(aes(x=log2FoldChange)) + geom_histogram()
  # 
  # FS12b@sam_data$AULC
  
  
  
  
  
  
  SCFA_OTU_assoc_treat <- 
    list(
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'butyrate', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'caproate', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'valerate', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'succinate', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'lactate', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'acetate', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'propionate', treatment = T)[[2]], 
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'isobutyrate', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'isovalerate', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'phenylacetate', treatment = T)[[2]],
      DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'oxalate', treatment = T)[[2]]
    ) %>% 
    bind_rows() %>%
    # filter(padj < 0.05 & abs(log2FoldChange) >.5) %>%
    mutate(node_name=paste(day, tissue, direction, covariate,sep = '_'))
    # 
  # bind_rows(SCFA_OTU_assoc_treat, sal_OTU_assoc_treat) %>%
  #   mutate(BH_all_P = p.adjust(pvalue, method = 'fdr')) %>%
  #   filter(abs(log2FoldChange) > .5) %>%
  #   ggplot(aes(x=BH_all_P, y=pvalue)) + geom_point()
  # 
  
  
  # # FDR on all tests for both scfa and sal
  edges <- bind_rows(SCFA_OTU_assoc_treat, sal_OTU_assoc_treat) %>% 
    mutate(BH_all_P = p.adjust(pvalue, method = 'fdr')) %>%   
    filter(BH_all_P < EDGE_QVAL) %>%
    filter(abs(log2FoldChange) > EDGE_LFC) %>%
    transmute(from=as.character(OTU), 
              to=sub('D[0-9]+_[A-Z]_','',node_name), 
              weight=abs(log2FoldChange))
  
  
  
  
  nodes <- tibble(
    V_ID=unique(c(edges$from, edges$to)), 
    type=case_when(
      grepl('Otu', V_ID) ~ 'OTU',
      grepl('ate$', V_ID) ~ 'SCFA',
      grepl('log_sal|AULC', V_ID) ~ 'Salmonella',
      TRUE ~ 'ERR'
    )
  )
  
  
  nodes[nodes$type == 'ERR',]
  
  NET <- fortify(as.edgedf(edges), nodes)
  
  ### geomnet
  set.seed(2)
  
  gg <- ggplot(data = NET, aes(from_id = from_id, to_id = to_id)) +
    geom_net(colour = "darkred", labelon=TRUE, size = 5,layout.alg = 'fruchtermanreingold', 
             directed = FALSE, vjust = 0.5,fontsize=2, labelcolour = "black",
             ecolour = "grey40") +
    theme_net()
  
  gg
  
  graph_layout <- 
    gg %>% ggplot_build()
  
  
  
  GRAPH_EDGES <- graph_layout$data[[1]] %>%  select(from , to, x, y, xend, yend)
  
  gn1TMP <- 
    GRAPH_EDGES %>% 
    select(from, x, y) %>%
    unique()
  
  gn2TMP <- 
    GRAPH_EDGES %>% 
    select(to, xend, yend) %>%
    unique() %>%
    transmute(from=to, x=xend, y=yend)
  
  NODE_DATA <- bind_rows(gn1TMP, gn2TMP) %>% unique()
  SCFA_NODES <- NODE_DATA %>% filter(!grepl('Otu', from)) # should be 'other_nodes'
  
  
  ONODES <- 
    NODE_DATA %>%
    filter(grepl('Otu', from)) %>% 
    transmute(OTU=from, x=x, y=y) %>% 
    left_join(taxtab)
  
  unique(ONODES$Family)
  
  # NOW BRING IN ABUND DATA  
  # modify this to be aggregate over all treatments
  
  ONODES <- 
    all_agg_abund_info %>% 
    # filter(treatment == 'RPS') %>% 
    filter(OTU %in%ONODES$OTU) %>% 
    select(OTU, perc_comm) %>% 
    right_join(ONODES) %>% 
    mutate(Class=fct_reorder(Class, perc_comm, .fun = sum, .desc = TRUE))
  
  ### MERGE IN TREATMENT ENRICHMENT HERE
  ONODES %>% dplyr::count(OTU) %>% filter(n!=1)
  
  unique(ONODES$Family)
  
  
  
  RPS_enrich_OTUs <- 
    tibble(OTU=RPS_sigOTUs, 
           enrich='RPS') %>% right_join(ONODES) %>% 
    filter(enrich == 'RPS')
  
  
  library(ggrepel)
  
  FIG_S7 <- 
    ONODES %>%
    ggplot(aes(x=x, y=y)) + 
    geom_point(data=RPS_enrich_OTUs, aes(x=x, y=y), color='blue', size=10, alpha=.5)+
    geom_segment(data=GRAPH_EDGES, aes(x=x, y=y, xend=xend, yend=yend), color='grey') + 
    geom_point(shape=21, aes(size=perc_comm, fill=Class)) + 
    scale_fill_brewer(palette = 'Set1') + 
    # geom_text(aes(label=Genus), size=3, nudge_y = .02) + 
    geom_point(data = SCFA_NODES, aes(x=x, y=y), size=3, shape=22, fill='black') + 
    geom_label(data = SCFA_NODES, aes(x=x, y=y, label=from), size=2,color='white', fill='black')  + 
    # theme_net() +
    geom_text_repel(aes(label=Genus),size=2, max.overlaps = 100000)+
    theme(legend.position = 'bottom') + 
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank())
  
  
  return(FIG_S7)
  
}

all_treats <- test_fun(phyloseq_obj = FS12b, EDGE_QVAL = .05, EDGE_LFC = .5)


just_RPS <- test_fun(phyloseq_obj = FS12_RPS, EDGE_QVAL = .05, EDGE_LFC = .5)

# all_melt <- 
#   FS12b %>% 
#   rarefy_even_depth() %>% 
#   psmelt()
# 
# taxtab <- as(tax_table(FS12b), 'matrix') %>% data.frame() %>% rownames_to_column(var = 'OTU')
# 
# all_agg_abund_info <- 
#   all_melt %>%
#   group_by(OTU) %>% 
#   summarise(mean_abund=mean(Abundance)) %>% 
#   mutate(tot_abund=sum(mean_abund), 
#          prop_comm=mean_abund/tot_abund, 
#          perc_comm=prop_comm * 100) %>% 
#   ungroup() %>% 
#   left_join(taxtab)
# 
# 
# all_agg_abund_info %>% dplyr::count(OTU) %>% filter(n !=1)
# 
# 
# # padj after all calcs?
# 
# sal_OTU_assoc_treat <- 
#   list(
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'log_sal', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'log_sal', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'X', covariate = 'log_sal', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D0', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D2', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D7', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D14', tissue = 'F', covariate = 'AULC', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'F', covariate = 'AULC', treatment = T)[[2]]
#     ) %>% 
#   bind_rows() %>%
#   # filter(padj < 0.05 & abs(log2FoldChange) >.5) %>%
#   # filter(padj < 0.05 ) %>% 
#   mutate(node_name=paste(day, tissue, direction, covariate,sep = '_'))
# 
# sal_OTU_assoc_treat %>% filter(covariate == 'AULC') %>% 
#   ggplot(aes(x=log2FoldChange)) + geom_histogram()
# 
# FS12b@sam_data$AULC
# 
# 
# 
# 
# 
# 
# SCFA_OTU_assoc_treat <- 
#   list(
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'butyrate', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'caproate', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'valerate', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'succinate', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'lactate', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'acetate', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'propionate', treatment = T)[[2]], 
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'isobutyrate', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'isovalerate', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'phenylacetate', treatment = T)[[2]],
#     DESeq_cov_asso(phyloseq_obj = FS12b, day = 'D21', tissue = 'C', covariate = 'oxalate', treatment = T)[[2]]
#   ) %>% 
#   bind_rows() %>%
#   # filter(padj < 0.05 & abs(log2FoldChange) >.5) %>%
#   mutate(node_name=paste(day, tissue, direction, covariate,sep = '_'))
# 
# bind_rows(SCFA_OTU_assoc_treat, sal_OTU_assoc_treat) %>%
#   mutate(BH_all_P = p.adjust(pvalue, method = 'fdr')) %>% 
#   filter(abs(log2FoldChange) > .5) %>% 
#   ggplot(aes(x=BH_all_P, y=pvalue)) + geom_point()
# 
# 
# #non OTU correlations 
# 
# 
# 
# FS12b@sam_data %>% as_tibble() %>% 
#   select(pignum, AULC, log_sal,day, tissue, treatment) %>%
#   filter(pignum == 181)
# 
# 
# 
# 
# edges <- bind_rows(SCFA_OTU_assoc_treat, sal_OTU_assoc_treat) %>% 
#   mutate(BH_all_P = p.adjust(pvalue, method = 'fdr')) %>%   # FDR on all tests for both scfa and sal
#   filter(BH_all_P < 0.05) %>% 
#   filter(abs(log2FoldChange) > .5) %>% 
#   transmute(from=as.character(OTU), 
#             to=sub('D[0-9]+_[A-Z]_','',node_name), 
#             weight=abs(log2FoldChange))
# 
# 
# 
# 
# nodes <- tibble(
#   V_ID=unique(c(edges$from, edges$to)), 
#   type=case_when(
#     grepl('Otu', V_ID) ~ 'OTU',
#     grepl('ate$', V_ID) ~ 'SCFA',
#     grepl('log_sal|AULC', V_ID) ~ 'Salmonella',
#     TRUE ~ 'ERR'
#   )
# )
# 
# 
# nodes[nodes$type == 'ERR',]
# 
# NET <- fortify(as.edgedf(edges), nodes)
# 
# ### geomnet
# set.seed(2)
# 
# gg <- ggplot(data = NET, aes(from_id = from_id, to_id = to_id)) +
#   geom_net(colour = "darkred", labelon=TRUE, size = 5,layout.alg = 'fruchtermanreingold', 
#            directed = FALSE, vjust = 0.5,fontsize=2, labelcolour = "black",
#            ecolour = "grey40") +
#   theme_net()
# 
# gg
# 
# graph_layout <- 
#   gg %>% ggplot_build()
# 
# 
# 
# GRAPH_EDGES <- graph_layout$data[[1]] %>%  select(from , to, x, y, xend, yend)
# 
# gn1TMP <- 
#   GRAPH_EDGES %>% 
#   select(from, x, y) %>%
#   unique()
# 
# gn2TMP <- 
#   GRAPH_EDGES %>% 
#   select(to, xend, yend) %>%
#   unique() %>%
#   transmute(from=to, x=xend, y=yend)
# 
# NODE_DATA <- bind_rows(gn1TMP, gn2TMP) %>% unique()
# SCFA_NODES <- NODE_DATA %>% filter(!grepl('Otu', from)) # should be 'other_nodes'
# 
# 
# ONODES <- 
#   NODE_DATA %>%
#   filter(grepl('Otu', from)) %>% 
#   transmute(OTU=from, x=x, y=y) %>% 
#   left_join(taxtab)
# 
# unique(ONODES$Family)
# 
# # NOW BRING IN ABUND DATA  
# # modify this to be aggregate over all treatments
# 
# ONODES <- 
#   all_agg_abund_info %>% 
#   # filter(treatment == 'RPS') %>% 
#   filter(OTU %in%ONODES$OTU) %>% 
#   select(OTU, perc_comm) %>% 
#   right_join(ONODES) %>% 
#   mutate(Class=fct_reorder(Class, perc_comm, .fun = sum, .desc = TRUE))
# 
# ### MERGE IN TREATMENT ENRICHMENT HERE
# ONODES %>% dplyr::count(OTU) %>% filter(n!=1)
# 
# unique(ONODES$Family)
# 
# 
# 
# RPS_enrich_OTUs <- 
#   tibble(OTU=RPS_sigOTUs, 
#          enrich='RPS') %>% right_join(ONODES) %>% 
#   filter(enrich == 'RPS')
# 
# 
# library(ggrepel)
# 
# FIG_S7 <- 
#   ONODES %>%
#   ggplot(aes(x=x, y=y)) + 
#   geom_point(data=RPS_enrich_OTUs, aes(x=x, y=y), color='blue', size=10, alpha=.5)+
#   geom_segment(data=GRAPH_EDGES, aes(x=x, y=y, xend=xend, yend=yend), color='grey') + 
#   geom_point(shape=21, aes(size=perc_comm, fill=Class)) + 
#   scale_fill_brewer(palette = 'Set1') + 
#   # geom_text(aes(label=Genus), size=3, nudge_y = .02) + 
#   geom_point(data = SCFA_NODES, aes(x=x, y=y), size=3, shape=22, fill='black') + 
#   geom_label(data = SCFA_NODES, aes(x=x, y=y, label=from), size=2,color='white', fill='black')  + 
#   # theme_net() + 
#   geom_text_repel(aes(label=Genus),size=2, max.overlaps = 100000)+
#   theme(legend.position = 'bottom')
# 
# 
# FIG_S7
# 
