library(ANCOMBC)


ancom_DA <- 
  ancombc(phyloseq = FS12b, formula = 'log_sal + treatment', p_adj_method = 'fdr')

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



ancom_subset <- 
  function(PHYLO, DAY, TISSUE, FORMULA, CONSERVE=F){
    # browser()
  keep_samps <- PHYLO@sam_data$tissue == TISSUE & PHYLO@sam_data$day == DAY
  filt_phylo <- prune_samples(keep_samps, PHYLO)
  
  ancom_DA <- 
    ancombc(phyloseq = filt_phylo,
            formula = FORMULA,
            p_adj_method = 'fdr', conserve = CONSERVE)
  
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



