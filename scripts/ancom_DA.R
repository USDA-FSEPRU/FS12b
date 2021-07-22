library(ANCOMBC)


ancom_DA <- 
  ancombc(phyloseq = FS12b,
          formula = 'log_sal + treatment', p_adj_method = 'fdr')

ancom_DA$res[['beta']] %>%
  rownames_to_column(var = 'OTU') %>%
  gather(key='covariate', value='value', -OTU)


tidy_param <- function(param_df, parameter_name){
  long_params <- 
    param_df %>%
    rownames_to_column(var='OTU') %>% 
    pivot_longer(values_to=parameter_name, names_to='covariate', -c('OTU'))
  return(long_params)
}


map2(.x=ancom_DA$res, .y=names(ancom_DA$res), tidy_param)%>%
  purrr::reduce(left_join) %>% mutate


#### functionize ###

FS12b@sam_data$treatment
FORMULA='log_sal + treatment'


FS12b@sam_data$day
FS12b@sam_data$tissue %>% unique()

ANCOM_DIFABUND_conserve <- 
  bind_rows(
  ancom_subset(PHYLO = FS12b, DAY = 'D0',TISSUE = 'F',FORMULA = 'treatment', CONSERVE = T),
  ancom_subset(PHYLO = FS12b, DAY = 'D2',TISSUE = 'F',FORMULA = 'treatment', CONSERVE = T),
  ancom_subset(PHYLO = FS12b, DAY = 'D7',TISSUE = 'F',FORMULA = 'treatment', CONSERVE = T),
  ancom_subset(PHYLO = FS12b, DAY = 'D14',TISSUE = 'F',FORMULA = 'treatment', CONSERVE = T),
  ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'F',FORMULA = 'treatment', CONSERVE = T), 
  ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'X',FORMULA = 'treatment', CONSERVE = T),
  ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'C',FORMULA = 'treatment', CONSERVE = T),
  ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'I',FORMULA = 'treatment', CONSERVE = T))



ANCOM_DIFABUND <- 
  bind_rows(
    ancom_subset(PHYLO = FS12b, DAY = 'D0',TISSUE = 'F',FORMULA = 'treatment'),
    ancom_subset(PHYLO = FS12b, DAY = 'D2',TISSUE = 'F',FORMULA = 'treatment'),
    ancom_subset(PHYLO = FS12b, DAY = 'D7',TISSUE = 'F',FORMULA = 'treatment'),
    ancom_subset(PHYLO = FS12b, DAY = 'D14',TISSUE = 'F',FORMULA = 'treatment'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'F',FORMULA = 'treatment'), 
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'X',FORMULA = 'treatment'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'C',FORMULA = 'treatment'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'I',FORMULA = 'treatment'))


ANCOM_DIFABUND_ZEROS <- 
  bind_rows(
    ancom_subset(PHYLO = FS12b, DAY = 'D0',TISSUE = 'F',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE),
    ancom_subset(PHYLO = FS12b, DAY = 'D2',TISSUE = 'F',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE),
    ancom_subset(PHYLO = FS12b, DAY = 'D7',TISSUE = 'F',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE),
    ancom_subset(PHYLO = FS12b, DAY = 'D14',TISSUE = 'F',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'F',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE), 
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'X',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'C',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'I',FORMULA = 'treatment', zero_group = 'treatment', STRUCT_ZERO = TRUE))


# ANCOM with zeros vs ANCOM w/o zeros

noZs <- ANCOM_DIFABUND %>% filter(q_val < 0.05) %>% write_tsv('./output/ANCOM_SIGs.tsv')

DESEQ_RES <- read_tsv('./output/Control_vs_All_DESeq.tsv')


DESEQ_RES <- 
  DESEQ_RES %>%
  transmute(OTU, DESeq_beta=log2FoldChange,
            day, tissue, treatment=sub('down_', '',Treatment)) %>% 
  transmute(OTU, DESeq_beta, day, tissue,
            DESeq_result=ifelse(DESeq_beta > 0 , 
                             paste('up', treatment),
                             paste('down', treatment)))
  



ANCOM_RES <- 
  noZs %>% 
  transmute(OTU, ANCOM_beta = beta, day, tissue, treatment=sub('treatment','',covariate)) %>% 
  transmute(OTU, ANCOM_beta, day, tissue, 
            ANCOM_result=ifelse(ANCOM_beta > 0 , 
                                paste('up', treatment),
                                paste('down', treatment)))

look <- ANCOM_RES %>%
  filter(grepl('RPS', ANCOM_result)) %>%
  full_join(DESEQ_RES %>% filter(grepl('RPS', DESeq_result))) %>% 
  mutate(DESeq_result=ifelse(is.na(DESeq_result), 'Not DA', DESeq_result)) %>% 
  mutate(AGREE=DESeq_result == ANCOM_result) 

look %>% filter(AGREE ==FALSE)

look %>% ggplot(aes(x=AGREE, y=ANCOM_beta, fill=DESeq_result)) + geom_point(shape=21)


look %>% ggplot(aes(x=AGREE, y=))

DESEQ_RES %>% 
  left_join(ANCOM_RES) %>%
  mutate(AGREE=DESeq_result == ANCOM_result)





#


Zs <- ANCOM_DIFABUND_ZEROS %>% filter(q_val < 0.05)

noZs %>% group_by(day, covariate, tissue) %>% tally() %>% arrange(desc(n))
Zs %>% group_by(day, covariate, tissue) %>% tally() %>% arrange(desc(n))
ancom_subset <- 
  function(PHYLO, DAY, TISSUE,
           FORMULA, CONSERVE=F,
           GLOBAL=F, zero_group=NULL, STRUCT_ZERO=FALSE){
    # browser()
  keep_samps <- PHYLO@sam_data$tissue == TISSUE & PHYLO@sam_data$day == DAY
  filt_phylo <- prune_samples(keep_samps, PHYLO)
  
  ancom_DA <- 
    ancombc(phyloseq = filt_phylo,
            formula = FORMULA,
            p_adj_method = 'fdr',
            conserve = CONSERVE, group = zero_group, struc_zero = STRUCT_ZERO,
            global=GLOBAL)
  
  map2(.x=ancom_DA$res, .y=names(ancom_DA$res), tidy_param)%>%
    purrr::reduce(left_join) %>%
    mutate(day=DAY, 
           tissue=TISSUE)
}


tidy_param <- function(param_df, parameter_name){
  long_params <- 
    param_df %>%
    rownames_to_column(var='OTU') %>% 
    pivot_longer(values_to=parameter_name, names_to='covariate', -c('OTU'))
  return(long_params)
}
# 
# 
# ANCOM_DIFABUND %>%
#   mutate(FDR_P=p.adjust(p_val, method = 'fdr')) %>%
#   arrange(FDR_P)
# 

taxtab_df <- 
  as(FS12b@tax_table, 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var='OTU')



ANCOM_DIFABUND_conserve %>%
  left_join(taxtab_df) %>% 
  filter(tissue =='F') %>%
  filter(q_val < 0.05) %>%
  filter(beta >0) %>% 
  filter(covariate == 'treatmentRPS') %>%
  # filter(day=='D0') %>% 
  ggplot(aes(x=Genus, y=beta, color=covariate, shape=day)) +
  geom_point() +
  coord_flip()
# 



ANCOM_DIFABUND_conserve %>%
  left_join(taxtab_df) %>% 
  filter(tissue =='F') %>%
  filter(q_val < 0.05) %>%
  filter(beta <0) %>% 
  filter(covariate == 'treatmentAcid') %>%
  # filter(day=='D0') %>% 
  ggplot(aes(x=Genus, y=beta, color=covariate, shape=day)) +
  geom_point() +
  coord_flip()
# 


#### log_sal assoc ###


ANCOM_log_sal <- 
  bind_rows(
    ancom_subset(PHYLO = FS12b, DAY = 'D2',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D7',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D14',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'F',FORMULA = 'treatment + log_sal'), 
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'X',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'C',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'I',FORMULA = 'treatment + log_sal'))


ANCOM_log_sal %>% filter(covariate == 'log_sal') %>% 
  filter(q_val < 0.05) %>% left_join(taxtab_df) %>% 
  ggplot(aes(x=Genus, y=beta, color=tissue, shape=day)) +
  geom_point() +
  coord_flip()



ANCOM_log_sal <- 
  bind_rows(
    ancom_subset(PHYLO = FS12b, DAY = 'D2',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D7',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D14',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'F',FORMULA = 'treatment + log_sal'), 
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'X',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'C',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'I',FORMULA = 'treatment + log_sal'))


ANCOM_log_sal %>% filter(covariate == 'log_sal') %>% 
  filter(q_val < 0.05) %>% left_join(taxtab_df) %>% 
  ggplot(aes(x=Genus, y=beta, color=tissue, shape=day)) +
  geom_point() +
  coord_flip()

FS12b@sam_data <- 
  FS12b@sam_data %>%
  as('data.frame') %>%
  mutate(across(ends_with('ate'), ~ scale(.x, center = F)),
         log_sal=scale(log_sal, center=F), 
         AULC=scale(AULC, center=F)) %>% sample_data()




FORM <- 'treatment + log_sal + AULC + butyrate + caproate + valerate + succinate + acetate + propionate'


ALL_SCFAS <- ancom_subset(PHYLO = FS12b, DAY='D21', TISSUE = 'C', FORMULA = FORM)


ALL_SCFAS %>%
  filter(q_val < 0.05) %>% 
  arrange(desc(beta)) %>% 
  filter(grepl('ate|AULC|log_sal', covariate)) %>% 
  group_by(covariate) %>% tally() %>% arrange(desc(n))

ANCOM_log_sal <- 
  bind_rows(
    ancom_subset(PHYLO = FS12b, DAY = 'D0',TISSUE = 'F',FORMULA = 'treatment + AULC'), 
    ancom_subset(PHYLO = FS12b, DAY = 'D2',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D7',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D14',TISSUE = 'F',FORMULA = 'treatment + log_sal'),
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'F',FORMULA = 'treatment + log_sal'), 
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'X',FORMULA = 'treatment + log_sal'), 
    ancom_subset(PHYLO = FS12b, DAY = 'D21',TISSUE = 'I',FORMULA = 'treatment + log_sal')
    )


ANCOM_log_sal %>%
  filter(tissue =='F') %>% 
  filter(q_val < 0.05) %>%
  filter(grepl('log_sal', covariate)) %>% 
  arrange(desc(beta)) %>%
  left_join(taxtab_df) %>% 
  ggplot(aes(x=Genus, y=beta, fill=day)) +
  geom_point(shape=21) + 
  facet_wrap(~covariate, scales = 'free', ncol=1) + 
  coord_flip()


FS12b_RPS <- FS12b %>% prune_samples(samples = FS12b@sam_data$treatment == 'RPS')

ANCOM_log_sal <- 
  bind_rows(
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D0',TISSUE = 'F',FORMULA = 'AULC'), 
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D2',TISSUE = 'F',FORMULA = 'log_sal'),
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D7',TISSUE = 'F',FORMULA = 'log_sal'),
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D14',TISSUE = 'F',FORMULA = 'log_sal'),
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D21',TISSUE = 'F',FORMULA = 'log_sal'), 
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D21',TISSUE = 'X',FORMULA = 'log_sal'), 
    ancom_subset(PHYLO = FS12b_RPS, DAY = 'D21',TISSUE = 'I',FORMULA = 'log_sal')
  ) 


ANCOM_log_sal <- ANCOM_log_sal %>%
  mutate(FDR_pval = p.adjust(p_val, method = 'fdr')) %>% 
  filter(FDR_pval < 0.05)


# Caclulate associations with cecal SCFAs
FORM <- 'caproate + valerate + butyrate + succinate + acetate + propionate'

ALL_SCFAS <- ancom_subset(PHYLO = FS12b_RPS, DAY='D21', TISSUE = 'C', FORMULA = FORM, CONSERVE = T)

ANCOM_SCFAS <- 
  ALL_SCFAS %>%
  mutate(FDR_pval = p.adjust(p_val, method = 'fdr')) %>% 
  filter(FDR_pval < 0.05)

ANCOM_SCFAS %>% ggplot(aes(x=FDR_pval, y=q_val)) + geom_point() + xlim(0,.1) + ylim(0,.1)

# ALL_SCFAS <- ALL_SCFAS %>% filter(q_val <0.05)
# FS12b_RPS %>% prune_samples(samples = FS12b_RPS@sam_data$day == 'D21') %>% sample_data()
# 
# 
# ANCOM_SCFAS %>%
#   filter(q_val < 0.05) %>% 
#   left_join(taxtab_df) %>% 
#   ggplot(aes(x=Genus, y=beta, fill=covariate)) + 
#   geom_point(shape=21) + 
#   coord_flip() + ylim(-1,1)
# 
# hist(ALL_SCFAS$beta, breaks = 100)

# 
# look <- 
#   ALL_SCFAS %>%
#   filter(q_val < 0.05) %>% 
#   left_join(taxtab_df) %>% 
#   filter(beta >0) %>% 
#   filter(covariate %in% c('butyrate', 'valerate', 'caproate'))


ANCOM_log_sal %>%
  filter(q_val < 0.05) %>%
  left_join(taxtab_df) %>% 
  ggplot(aes(x=Genus, y=beta, fill=day)) +
  geom_point(shape=21) + 
  coord_flip()

sal_OTUS <- ANCOM_log_sal %>% pull(OTU) %>% unique()


SCFA_OTUS <- ANCOM_SCFAS %>% pull(OTU) %>% unique() 

sum(sal_OTUS %in% SCFA_OTUS)

keepers <- sal_OTUS[sal_OTUS %in% SCFA_OTUS]



EDGE_INFO <- 
  bind_rows(ANCOM_log_sal, ANCOM_SCFAS) %>% 
  filter(OTU %in% keepers) %>% 
  filter(tissue != 'I') %>%
  mutate(NODE_NAME=ifelse(beta >0,
                          paste('increased', covariate),
                          paste('decreased', covariate)), 
         edge_weight=abs(beta)) %>% 
  left_join(taxtab_df)



edge_df <- 
  EDGE_INFO %>% 
  transmute(from=OTU, to=NODE_NAME, edge_weight = edge_weight, edge_pval=FDR_pval)

node_df <- 
  bind_rows(EDGE_INFO %>% transmute(NODE_NAME = OTU, 
                                    type='OTU'), 
            EDGE_INFO %>% transmute(NODE_NAME = NODE_NAME,
                                    type=ifelse(grepl('ate$',covariate), 'SCFA', 'Salmonella')))





node_df %>% pull(type) %>% unique()


library(geomnet)


