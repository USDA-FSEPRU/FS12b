
library(tidyverse)
library(phyloseq)
library(vegan)
library(funfuns)
library(broom)
library(DESeq2)
library(cowplot)
library(geomnet)
library(ggrepel)


DESeq_cov_asso <- 
  function(phyloseq_obj,
           day,
           tissue,
           covariate,
           shrink_type='apeglm',
           cookscut=TRUE, 
           treatment=FALSE, 
           scale_cov=TRUE, 
           p.plot=0.05, 
           plot_lab='Genus', 
           l2fc_plot=0.25){
    
    num_treats <- length(unique(phyloseq_obj@sam_data[['treatment']]))
    # include treatment in formula?
    if (treatment & (num_treats > 1)){
      
      form <- formula(paste('~', covariate, '+ treatment'))
      print(form)
    } else {
      form <- formula(paste('~', covariate))
      print(form)
    }
    
    
    # print(form)
    # browser()
    FS12b.glom <- phyloseq_obj %>%
      prune_samples(samples = phyloseq_obj@sam_data$day %in% c(day) & phyloseq_obj@sam_data$tissue == tissue & !is.na(phyloseq_obj@sam_data[[covariate]]))
    
    # scale covariate?
    if (scale_cov){
      print('scale_cov == TRUE, scaling covariate')
      FS12b.glom@sam_data[[covariate]] <- scale(FS12b.glom@sam_data[[covariate]])
      
    }
    
    
    FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
    
    # FS12b.glom@sam_data$log_sal
    
    FS12b.de <- phyloseq_to_deseq2(FS12b.glom, form)
    FS12b.de <- DESeq(FS12b.de, test = 'Wald', fitType = 'parametric')
    
    # these are not both possible.  Right now only lfcshrink is doing anytihng
    # res <- results(FS12b.de, cooksCutoff = FALSE, name = covariate)
    res <- results(FS12b.de, name=covariate, cooksCutoff = cookscut)
    sigtab <- lfcShrink(dds = FS12b.de, res=res, coef = covariate, type = shrink_type)
    
    # browser()
    # resultsNames(FS12b.de)
    
    
    # res <- res[!is.na(res$padj),]
    # res <- res[res$padj < 0.1,]
    # sigtab <- res[abs(res$log2FoldChange) > .1 ,]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
    sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = TRUE)
    # sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
    sigtab$OTU <- rownames(sigtab)
    sigtab[['direction']] <- ifelse(sigtab$log2FoldChange > 0 , 'increased', 'decreased')
    # sigtab$salm <- ifelse(sigtab$log2FoldChange >0 , 'increased', 'decreased')
    sigtab <- sigtab[order(sigtab$log2FoldChange),]
    sigtab$OTU <- factor(sigtab$OTU, levels = sigtab$OTU)
    sigtab$day <- day
    sigtab$tissue <- tissue
    sigtab[['covariate']] <- covariate
    plot_tab <- sigtab %>% filter(padj < p.plot & abs(log2FoldChange) > l2fc_plot)
    fig_tit <- paste(covariate, 'associations with OTUs')
    if (nrow(plot_tab) > 0 ){
      p <- plot_tab %>% 
        ggplot(aes_string(x='OTU', y='log2FoldChange', fill='direction')) +
        geom_col(color='black') +
        coord_flip() +
        geom_text(aes_string(label=plot_lab, y=0)) + 
        ggtitle(fig_tit)
      
    } else {
      print('nothing to plot')
      p <- ggplot() + 
        ggtitle('No sig associations')
    }
    
    return(list(p, sigtab))
    
    
  }
# source('./scripts/FS12b_16s_figs.R')


RPS_enriched_OTUs <- 
  read_tsv('./output/Control_vs_All_DESeq.tsv') %>%
  filter(Treatment== 'RPS') %>% 
  pull(OTU) %>% 
  unique()


## Build phyloseq object ##
## probably split this into the calculations and then the figures.
### write out the calcs and then read back in in figures script.


# OTU <- 
#   read_tsv('./data/FS12b_OTU.tsv') %>%
#   column_to_rownames(var='Group') %>%
#   as.matrix() %>% 
#   otu_table(taxa_are_rows = FALSE)
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
    mutate(Class=fct_reorder(Class, perc_comm, .fun = sum, .desc = TRUE)) %>% 
    mutate(Class=fct_lump_n(Class, 9, w=perc_comm))
  
  ### MERGE IN TREATMENT ENRICHMENT HERE
  # ONODES %>% dplyr::count(OTU) %>% filter(n!=1)
  
  # unique(ONODES$Family)
  
  
  
  RPS_enrich_OTUs <- 
    tibble(OTU=highlight_these_otus, 
           enrich='RPS') %>%
    right_join(ONODES) %>% 
    filter(enrich == 'RPS')
  
  
  
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

FS12b_RPS <- FS12b %>% prune_samples(samples=FS12b@sam_data$treatment == 'RPS')


all_treats <- test_fun(phyloseq_obj = FS12b, 
                       EDGE_QVAL = .05,
                       EDGE_LFC = .5,
                       SEED = 2,
                       highlight_these_otus = RPS_enriched_OTUs)


just_RPS <- test_fun(phyloseq_obj = FS12b_RPS, 
                     EDGE_QVAL = .05,
                     EDGE_LFC = .5, 
                     SEED = 2, 
                     highlight_these_otus = RPS_enriched_OTUs)

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




build_net <- 
  function(phyloseq_obj, EDGE_QVAL, EDGE_LFC, SEED, highlight_these_otus){
    # browser()
    
    set.seed(SEED)
    FS12b <- phyloseq_obj
    
    
    
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
    return(list(NET, edges, nodes))
  }



plt_net <- 
  function(phyloseq_obj,
           NODES, EDGES, 
           highlight_these_nodes,
           LAYOUT='kamadakawai', 
           layout.par=list(NULL), 
           SEED=2){
    # NET needs to be a fortified edgelist
  # wrap this
  # plot network
    # browser()
    set.seed(SEED)
  all_melt <- 
    phyloseq_obj %>% 
    rarefy_even_depth() %>% 
    psmelt()
  
  taxtab <- as(tax_table(phyloseq_obj), 'matrix') %>%
    data.frame() %>%
    rownames_to_column(var = 'OTU')
  
  all_agg_abund_info <- 
    all_melt %>%
    group_by(OTU) %>% 
    summarise(mean_abund=mean(Abundance)) %>% 
    mutate(tot_abund=sum(mean_abund), 
           prop_comm=mean_abund/tot_abund, 
           perc_comm=prop_comm * 100) %>% 
    ungroup() %>% 
    left_join(taxtab)
  
  # removes edges and nodes relating to OTUs not in the abund info dfs
  EDGES <- EDGES[EDGES$from %in% all_agg_abund_info$OTU,]
  NODES <- NODES[grepl('creased', NODES$V_ID) | NODES$V_ID %in% all_agg_abund_info$OTU,]
  
 NET <-  fortify(as.edgedf(EDGES), NODES)
  ### geomnet
  
  
  gg <- ggplot(data = NET, aes(from_id = from_id, to_id = to_id)) +
    geom_net(colour = "darkred", labelon=TRUE, size = 5,layout.alg = LAYOUT, layout.par=layout.par,
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
    # filter(!is.na(perc_comm)) %>% 
    mutate(Class=fct_reorder(Class, perc_comm, .fun = sum, .desc = TRUE)) %>% 
    mutate(Class=fct_lump_n(Class, 8, w=perc_comm))
  
  ### MERGE IN TREATMENT ENRICHMENT HERE
  # ONODES %>% dplyr::count(OTU) %>% filter(n!=1)
  
  # unique(ONODES$Family)
  
  
  
  RPS_enrich_OTUs <- 
    tibble(OTU=highlight_these_nodes, 
           enrich='RPS') %>%
    right_join(ONODES) %>% 
    filter(enrich == 'RPS')
  
  
  
  FIG_S7 <- 
    ONODES %>%
    ggplot(aes(x=x, y=y)) + 
    geom_point(data=RPS_enrich_OTUs, aes(x=x, y=y), color='blue', size=10, alpha=.5)+
    geom_segment(data=GRAPH_EDGES, aes(x=x, y=y, xend=xend, yend=yend), color='grey') + 
    geom_point(shape=21, aes(size=perc_comm, fill=Class)) + 
    scale_fill_brewer(palette = 'Set1') + 
    # geom_text(aes(label=Genus), size=3, nudge_y = .02) + 
    # geom_point(data = SCFA_NODES, aes(x=x, y=y), size=3, shape=22, fill='black') + 
    geom_label(data = SCFA_NODES, aes(x=x, y=y, label=from), size=2,color='white', fill='black')  + 
    # theme_net() +
    geom_text_repel(aes(label=Genus),size=2.5, max.overlaps = 100000)+
    theme(legend.position = 'bottom') + 
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank())
  
  
  return(FIG_S7)
  
  }



#### All treatments together ##

full_net <- build_net(phyloseq_obj = FS12b,
                          EDGE_QVAL = 0.05, 
                          EDGE_LFC = .5, 
                          SEED = 2,
                          highlight_these_otus = RPS_enriched_OTUs)




### FILTER SMALL UNCONNECTED COMMUNITIES

g <- igraph::graph_from_data_frame(vertices = full_net[[3]], d=full_net[[2]], directed = F)
set.seed(2)
clouv <- cluster_louvain(g)

CONNECTED <- membership(clusters(g))


keep_these_nodes <- names(CONNECTED[CONNECTED %in% c(1, 2)])



filt_net <- 
  full_net[[1]] %>% 
  filter(from_id %in% keep_these_nodes) %>% 
  filter(to_id %in% keep_these_nodes)

filt_nodes <- 
  full_net[[3]] %>% filter(V_ID %in% keep_these_nodes)

filt_edges <- full_net[[2]] %>%
  filter(from %in% keep_these_nodes) %>%
  filter(to %in% keep_these_nodes)


##


full_net_fig7 <- 
  plt_net(phyloseq_obj = FS12b,
          NODES = filt_nodes, 
          EDGES = filt_edges, 
          highlight_these_nodes = RPS_enriched_OTUs,
          LAYOUT ='kamadakawai', layout.par = list(niter=10000))
full_net_fig7


ggsave('./output/full_net_fig7.jpeg', width = 9, height = 7, units = 'in')





### connections in just RPS treatment #



RPS_full_net <- build_net(phyloseq_obj = FS12b_RPS,
                EDGE_QVAL = 0.05, 
                EDGE_LFC = .5, 
                SEED = 2,
                highlight_these_otus = RPS_enriched_OTUs)


### FILTER THE SMALL ONES FROM THIS ###
# RPS_full_net_supp_fig <- plt_net(phyloseq_obj = FS12b,
#                              NODES = RPS_full_net[[3]],
#                              EDGES = RPS_full_net[[2]], 
#                              highlight_these_nodes = RPS_enriched_OTUs,
#                              LAYOUT ='kamadakawai', 
#                              layout.par = list(niter=10000))

library(igraph)


g <- igraph::graph_from_data_frame(vertices = RPS_full_net[[3]], d=RPS_full_net[[2]], directed = F)
set.seed(2)
clouv <- cluster_louvain(g)

CONNECTED <- membership(clusters(g))


keep_these_nodes <- names(CONNECTED[CONNECTED == 1])



filt_net <- 
  RPS_full_net[[1]] %>% 
  filter(from_id %in% keep_these_nodes) %>% 
  filter(to_id %in% keep_these_nodes)

filt_nodes <- 
  RPS_full_net[[3]] %>% filter(V_ID %in% keep_these_nodes)

filt_edges <- RPS_full_net[[2]] %>%
  filter(from %in% keep_these_nodes) %>%
  filter(to %in% keep_these_nodes)

RPS_full_net_supp_fig <- 
  plt_net(phyloseq_obj = FS12b_RPS,
          NODES = filt_nodes, 
          EDGES = filt_edges, 
          highlight_these_nodes = RPS_enriched_OTUs,
          LAYOUT ='kamadakawai', layout.par = list(niter=10000))
RPS_full_net_supp_fig


ggsave('./output/RPS_full_net_supp_fig.jpeg', width = 9, height = 7, units = 'in')



#####
# for fig 8
mem <- membership(clouv)

plot(clouv, g)

hist(degree(g))

DEGREE <- degree(g)



hist(DEGREE[grepl('Otu', names(DEGREE))])
otu_degrees <- DEGREE[grepl('Otu', names(DEGREE))]




central_OTUs <- otu_degrees[otu_degrees > 3] %>% sort(decreasing = T)
RPS_tax_table <- FS12b_RPS@tax_table %>%
  as.data.frame() %>%
  rownames_to_column(var='OTU')

tibble(OTU=names(central_OTUs), 
       DEGREE=central_OTUs) %>% 
  left_join(RPS_tax_table)




keeps <- names(mem) %in% c('increased_butyrate',
                  'increased_valerate',
                  'increased_caproate',
                  'decreased_log_sal', 
                  'increased_succinate', 
                  'decreased_AULC')

mem[keeps] %>% unique()

# 5, 7, 1

keepers <- membership(clouv)[(membership(clouv) %in% c(1,5,7))] %>% names()

filt_net <- 
  RPS_full_net[[1]] %>% 
  filter(from_id %in% keepers) %>% 
  filter(to_id %in% keepers)

filt_nodes <- 
  RPS_full_net[[3]] %>% filter(V_ID %in% keepers)

filt_edges <- RPS_full_net[[2]] %>%
  filter(from %in% keepers) %>%
  filter(to %in% keepers)

fig8 <- plt_net(phyloseq_obj = FS12b_RPS,
               NODES = filt_nodes, 
               EDGES = filt_edges, 
               highlight_these_nodes = RPS_enriched_OTUs,
               LAYOUT ='kamadakawai', layout.par = list(niter=10000))
fig8


ggsave('./output/figure8.jpeg', width = 9, height = 7, units = 'in')
# 
# plt_net(phyloseq_obj = FS12b_RPS,
#         NODES = BN[[3]],
#         EDGES = BN[[2]], 
#         highlight_these_nodes = RPS_enriched_OTUs,
#         LAYOUT = 'kamadakawai',
#         layout.par = list(niter=10000))
# 
# 
# 
# 
# 
# tstF <- plt_net(phyloseq_obj = FS12b,
#                NODES = filt_nodes, 
#                EDGES = filt_edges, 
#                highlight_these_nodes = RPS_enriched_OTUs,
#                LAYOUT ='fruchtermanreingold', layout.par = list(niter=10000))
# 
# tstF

