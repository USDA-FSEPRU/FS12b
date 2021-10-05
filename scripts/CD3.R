library(tidyverse)


CD3 <- read_csv('data/FS12b_CD3_Freq.csv')

CD3_long <- CD3 %>% 
  pivot_longer(cols = -c(`Flow ID`:Trt)) %>% 
  filter(grepl('CD4.CD8a.Foxp3.CD25.',name))

CD3_long %>%
  ggplot(aes(x=Trt, y=value, fill=Trt)) + 
  geom_boxplot() +
  facet_wrap(~name, scales = 'free')



CD3_long %>% group_by(Pig) %>% 
  summarise(tot=sum(value))


CD3_dat <- 
  CD3_long %>%
  mutate(pignum=Pig) %>% 
  select(-`Flow ID`, -Trt, -Pig) %>% 
  spread(key = name, value = value) 

# 
# CD3_corrs <- Hmisc::rcorr(CD3_mat) %>% 
#   broom::tidy() %>%
#   filter(p.value < 0.05) %>% 
#   filter(estimate > 0)









##### Vfas #######

# vfas <- read_tsv('./data/FS12b_meta.tsv')


colnames(vfas)
vfas <- read_csv('./data/FS12b_vfas.csv')
vfas <- vfas %>% filter(treatment %in% c('control', 'RPS', 'Acid','RCS'))
vfas$treatment <- factor(vfas$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

vfas.gather <- vfas %>% gather(key = VFA, value = mM, -(pignum:hour))

vfas.gather$set <- paste(vfas.gather$hour, vfas.gather$treatment)

scfas <- vfas %>% filter(hour == 0)



MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  mutate(treatment=factor(treatment, levels = c('Control', 'RPS', 'Acid', 'RCS')))


CD3_SCFA_SAL <- 
  MET %>%
  select(pignum,  AULC, log_sal, tissue, day) %>% 
  filter(day =='D21' & tissue == 'C') %>%
  right_join(scfas) %>% 
  left_join(CD3_dat) %>% 
  # filter(treatment == 'RPS') %>% 
  select(-tissue, -day, -treatment, -hour) %>% 
  column_to_rownames('pignum') %>% 
  as.matrix()



CD3_SCFA_SAL_COR <- 
  CD3_SCFA_SAL %>%
  Hmisc::rcorr() %>%
  broom::tidy() %>%
  filter(p.value < 0.05)


 





  # filter(treatment == 'RPS') %>% 
  # select(-tissue, -day, -hour, -treatment) %>%
  # column_to_rownames(var = 'pignum') %>% 
  # as.matrix()





