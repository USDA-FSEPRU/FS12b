# setwd('')
library(tidyverse)
library(cowplot)
library(broom)
theme_set(theme_cowplot())

vfas <- read_tsv('./data/FS12b_meta.tsv')


colnames(vfas)
vfas <- read_csv('./data/FS12b_vfas.csv')
vfas <- vfas %>% filter(treatment %in% c('Control', 'RPS', 'Acid','RCS'))
vfas$treatment <- factor(vfas$treatment, levels = c('Control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

vfas.gather <- vfas %>% gather(key = VFA, value = mM, -(pignum:hour))

vfas.gather$set <- paste(vfas.gather$hour, vfas.gather$treatment)

# initial peek #
# vfas %>% filter(pignum %in% c(458, 461, 469, 472))

### VFA CECUM TESTS ###

scfas <- vfas %>% filter(hour == 0)

but_test <- aov(butyrate~treatment, data = scfas)
cap_test <- aov(caproate~treatment, data = scfas)
val_test <- aov(valerate~treatment, data = scfas)
suc_test <- aov(succinate~treatment, data = scfas)

library(emmeans)


scfa_means <- 
  bind_rows(list(
  emmeans(but_test, specs =  ~ treatment) %>% tidy(conf.int=TRUE) %>% mutate(scfa='butyrate'),
  emmeans(val_test, specs =  ~ treatment) %>% tidy(conf.int=TRUE) %>% mutate(scfa='valerate'),
  emmeans(cap_test, specs =  ~ treatment) %>% tidy(conf.int=TRUE) %>% mutate(scfa='caproate'),
  emmeans(suc_test, specs =  ~ treatment) %>% tidy(conf.int=TRUE) %>% mutate(scfa='succinate'))) %>% 
  mutate(treatment=factor(treatment, levels = c('Control', 'RPS', 'Acid', 'RCS')), 
         scfa=factor(scfa, levels = c('butyrate', 'valerate', 'caproate', 'succinate')))

# summary(but_test)
# summary(val_test)
# summary(cap_test)
# summary(tot_test)

scfa_tests <- 
  bind_rows(list(
  TukeyHSD(but_test) %>% tidy() %>% mutate(scfa='butyrate'),
TukeyHSD(val_test) %>% tidy()%>% mutate(scfa='valerate'),
TukeyHSD(cap_test) %>% tidy()%>% mutate(scfa='caproate'),
TukeyHSD(suc_test) %>% tidy()%>% mutate(scfa='succinate')
)) %>% 
  filter(contrast %in% c('RCS-Control', 
                         'Acid-Control', 
                         'RPS-Control')) %>% 
  mutate(contrast=factor(contrast, levels=c('RCS-Control', 
                                            'Acid-Control', 
                                            'RPS-Control')),
         scfa=factor(scfa, levels = c('butyrate', 'caproate', 'valerate', 'succinate'))) %>% 
  mutate(p.plot=ifelse(adj.p.value < 0.05, adj.p.value, NA))


FIG6A <- 
  scfa_means %>%
  ggplot(aes(x=treatment, y=estimate, ymin=conf.low, ymax=conf.high))+
  geom_col(aes(fill=treatment))+
  geom_errorbar(width=.2) + 
  facet_wrap(~scfa, scales = 'free', nrow=1)+
  theme(legend.position = 'top') + 
  ylab('Cecal Concentation (mM)') + 
  xlab('') + 
  theme(axis.text.x=element_text(size=10, angle = -40, hjust=.1))+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))

FIG6A

FIG6B <- 
  scfa_tests %>%
  mutate(scfa=factor(scfa, levels = c('butyrate', 'valerate', 'caproate', 'succinate'))) %>% 
  filter(grepl('Control', contrast)) %>% 
  ggplot(aes(x=contrast, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(color=contrast)) + 
  geom_text(aes(label=round(p.plot, 3), color=contrast), nudge_x = .2)+
  coord_flip()+
  facet_wrap(~scfa, scales = 'free_x', nrow = 1) +
  theme(legend.position = 'none', 
        axis.text = element_text(size=11), 
        axis.title.x=element_text(size=11)) + 
  xlab('') + 
  ylab('difference from controls (mM)')+
  scale_color_manual(values=c('red', 'orange','#3399FF', 'red', 'grey', 'purple'))

FIG6B

ggdraw() + 
  draw_plot(FIG6A, x=0, y=.5, width = 1, height = .5) + 
  draw_plot(FIG6B, x=0, y=0, width = 1, height=.5)


## AULC VFA CORRELATIONS ##



MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  mutate(treatment=factor(treatment, levels = c('Control', 'RPS', 'Acid', 'RCS')))


RPS_cec_scfas <- 
  MET %>%
  select(pignum,  AULC, log_sal, tissue, day) %>% 
  filter(day =='D21' & tissue == 'C') %>%
  right_join(scfas) %>% 
  filter(treatment == 'RPS') %>%
  select(-tissue, -day, -hour, -treatment) %>%
  column_to_rownames(var = 'pignum') %>% 
  as.matrix()



cor_tests <- 
  Hmisc::rcorr(RPS_cec_scfas) %>%
  tidy() %>% 
  filter(column1 %in% c('AULC', 'butyrate', 'valerate', 'caproate', 'succinate')) %>% 
  filter(column2 %in% c('AULC', 'butyrate', 'valerate', 'caproate', 'succinate')) %>% 
  mutate(scfa=column1) %>% 
  filter(column2=='AULC') %>% 
  mutate(concentration=c(17, 5, .7, .45), 
         AULC=c(10,10,10,10), 
         stats=paste('r=', round(estimate, 2), '\nP=', round(p.value, 2)), 
         scfa=factor(scfa, levels = c('butyrate', 'valerate', 'caproate', 'succinate')))
  

FIG6C <- 
RPS_cec_scfas %>%
  as.data.frame() %>% 
  gather(-AULC, key='scfa', value='concentration') %>% 
  filter(scfa %in% c('butyrate', 'valerate', 'caproate', 'succinate')) %>%
  mutate(scfa=factor(scfa, levels = c('butyrate', 'valerate', 'caproate', 'succinate'))) %>% 
  ggplot(aes(x=concentration, y=AULC)) +
  geom_smooth(method = 'lm',se=F, color='black') + 
  geom_text(data=cor_tests, aes(label=stats), size=4, hjust=-.5)+
  geom_point(shape=21, size=3, fill='#3399FF', color='black') +
  facet_wrap(~scfa, scales = 'free')



fig6 <- ggdraw()+
  draw_plot(FIG6A, x = 0, y =.4, width = .6, height = .6)+
  draw_plot(FIG6B, x = 0, y = 0, width = .6, height = .45)+
  draw_plot(FIG6C, x = .6, y = 0, width = .4, height = 1)

fig6

ggsave(filename = './output/figure6.jpeg', device = 'jpeg', 
       width=280, height=180, units = 'mm')

### D21 CECUM ALL VFAS  ###




RPS_cec_scfas <- 
  MET %>%
  select(pignum,  AULC, log_sal, tissue, day) %>% 
  filter(day =='D21' & tissue == 'C') %>%
  right_join(scfas) %>% 
  # filter(treatment == 'RPS') %>%
  select(-tissue, -day, -hour, -treatment) %>%
  column_to_rownames(var = 'pignum') %>% 
  as.matrix()


RPS_cec_scfas_DF <-   MET %>%
  select(pignum,treatment,  AULC, log_sal, tissue, day) %>% 
  filter(day =='D21' & tissue == 'C') %>%
  right_join(scfas) %>% 
  # filter(treatment == 'RPS') %>%
  select(-tissue, -day, -hour)


cor_tests <- 
  Hmisc::rcorr(RPS_cec_scfas) %>%
  tidy() %>% 
  filter(column1 %in% c('AULC', 'butyrate', 'valerate', 'caproate', 'succinate')) %>% 
  filter(column2 %in% c('AULC', 'butyrate', 'valerate', 'caproate', 'succinate')) %>% 
  mutate(scfa=column1) %>% 
  filter(column2=='AULC') %>% 
  mutate(concentration=c(17, 5, .7, .45), 
         AULC=c(10,10,10,10), 
         stats=paste('r=', round(estimate, 2), '\nP=', round(p.value, 3)), 
         scfa=factor(scfa, levels = c('butyrate', 'valerate', 'caproate', 'succinate')))


TST <- 
  RPS_cec_scfas_DF %>%
  as.data.frame() %>% 
  gather(-c(AULC, pignum, treatment), key='scfa', value='concentration') %>% 
  filter(scfa %in% c('butyrate', 'valerate', 'caproate', 'succinate')) %>%
  mutate(scfa=factor(scfa, levels = c('butyrate', 'valerate', 'caproate', 'succinate'))) %>% 
  ggplot(aes(x=concentration, y=AULC)) +
  geom_smooth(method = 'lm',se=F, color='black') + 
  geom_text(data=cor_tests, aes(label=stats), size=4, hjust=-.5)+
  geom_point(aes(fill=treatment), shape=21, size=3, color='black') +
  facet_wrap(~scfa, scales = 'free')+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  
 

TST 

ggsave('./output/SCFA_AULC_cor_supp.jpg', height = 5, width = 7, units = 'in')
# #########
# 
# filter(vfas.gather, hour == 0 ) %>% 
#   ggplot(aes(x=treatment, y=mM, group=set, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +geom_jitter(shape=21, size=1.5, stroke=1, alpha=.75, width = .25)+
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   facet_wrap(~VFA, scales = 'free')
# 
# 
# unique(vfas.gather$VFA)
# 
# filter(vfas.gather, hour == 0 & VFA %in% c('butyrate', 'valerate', 'caproate', 'total')) %>%
#   ggplot(aes(x=treatment, y=mM, group=set, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +geom_jitter(shape=21, size=1.5, stroke=1, alpha=.75, width = .25)+
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   facet_wrap(~VFA, scales = 'free', nrow = 1) + ggtitle("Cecal SCFAs 21 days post challenge, no incubation")

#
# filter(vfas.gather, hour == 0 &
#          VFA %in% c('acetate', 'butyrate', 'propionate', 'valerate', 'caproate', 'total') & 
#          treatment %in% c('control', 'RPS')) %>%
#   ggplot(aes(x=treatment, y=mM, group=set, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +geom_jitter(shape=21, size=1.5, stroke=1, alpha=.75, width = .25)+
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   facet_wrap(~VFA, scales = 'free') + ggtitle("Cecal SCFAs 21 days post challenge")




# filter(vfas.gather, hour == 24) %>% ggplot(aes(x=treatment, y=mM, group=set, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +geom_jitter(shape=21, size=1.5, stroke=1, alpha=.75, width = .25)+
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   facet_wrap(~VFA, scales = 'free') + ggtitle("Cecal SCFAs 21 days post challenge, with incubation")

### propionate in 24hr incubation seems to have suffered in ZN

# 0 and 24 hr together

# filter(vfas.gather, VFA != 'succinate') %>% ggplot(aes(x=hour, y=mM, group=set, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   facet_wrap(~VFA, scales = 'free')
# 
###

# tests for RPS vs control

# buttest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'butyrate')
# valtest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'valerate')
# captest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'caproate')
# tottest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'total')
# 
# wilcox.test(buttest$mM~buttest$treatment) # p-value = 0.00399
# wilcox.test(valtest$mM~valtest$treatment) # p-value = 0.0007816
# wilcox.test(captest$mM~captest$treatment) # p-value = 0.01061
# wilcox.test(tottest$mM~tottest$treatment) # p-value = 0.1135
# 

### for correlations ###

# vfas0 <- filter(vfas, hour == 0) %>% select(-hour)
# 
# vfas24 <- filter(vfas, hour == 24) %>% select(-hour)
# 
# colnames(vfas24) <- paste('P', colnames(vfas24), sep = '_') # why 'P' ?  Because it shows the 'P'otential amount of VFAs for that bolus
# 
# colnames(vfas24)[1] <- 'pignum'
# 
# vfas_for_cor <- merge(vfas0, vfas24, by = 'pignum')
# 
# vfas_for_cor$day <- 21
# vfas_for_cor$tissue <- 'cecal_contents'
# vfas_for_cor <- vfas_for_cor %>% select(pignum, treatment, day, tissue, everything())


# write_csv(vfas_for_cor, './data/FS12b_vfas_for_cor.csv')

############### FECAL DATA #################




# library(tidyverse)


####### WARNING! I THINK THERE IS A PROBLEM WITH THE DATA FROM DAY 7 and 21!!!!!

### SUSPICIOUS AMOUNTS OF LACTATE AND SUCCINATE ARE PRESENT, USUALLY THESE ARE VERY RARE
### I THINK THERE WAS AN ISSUE DURING SAMPLE PREP
### MY GUESS IS THAT THE SAMPLES WERE ALLOWED TO THAW AT ROOM TEMP FOR TOO LONG AND 
### SOME BACTERIAL METABOLIC ACTIVITY GENERATED EXCESS LACTATE AND SUCCINATE
### USUALLY THESE WOULD BE CONVERTED TO OTHER SCFAS BY ANAEROBIC FERMENTATION BUT
### THE ANAEROBES PROBABLY DIED DUE TO OXYGEN EXPOSURE.


# write_csv(res.all, 'FS12b_fecal_vfas.csv')

# res.all <- read_csv('./data/FS12b_fecal_vfas.csv')
# # tests <- filter(res.all, time == 0) %>%
# #   do(pwilx=pairwise.wilcox.test(.$butyrate, .$treatment, p.adjust.method = 'none'))
# # #pairwise.wilcox.test(p.adjust.method = )
# # tests[[1]]
# # 
# # 
# # tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$caproate, .$treatment, p.adjust.method = 'none'))
# # #pairwise.wilcox.test(p.adjust.method = )
# # tests[[1]]
# # 
# # 
# # tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$valerate, .$treatment, p.adjust.method = 'none'))
# # #pairwise.wilcox.test(p.adjust.method = )
# # tests[[1]]
# # 
# # 
# # res.all$treatment <- factor(res.all$treatment, levels=c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))
# # 
# # 
# # res.all %>% ggplot(aes(x=treatment, y=acetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot() + geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=propionate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=butyrate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=valerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=caproate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=isovalerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=oxalate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=phenylacetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=succinate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=fumarate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # res.all %>% ggplot(aes(x=treatment, y=lactate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
# # 
# # 
# # res.gather <- res.all %>% gather(key = variable, value = concentration, -(pignum:set))
# # 
# # # res.gather %>% filter(variable == 'total')%>% ggplot(aes(x=time, y=concentration, group=pignum, color=treatment)) + geom_line()
# # 
# # res.sum <- res.gather %>%  group_by(time, treatment, variable) %>%
# #   summarise(mean=mean(concentration), sd=sd(concentration), n=n(), sterr=sd/sqrt(n))
# # 
# # res.sum <- res.sum %>% mutate(treat_var=paste(treatment, variable, sep='_'))
# # 
# # res.sum$treatment <- factor(res.sum$treatment, levels= c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))
# # 
# # res.sum %>% ggplot(aes(x=time, y=mean, group=treat_var, color=treatment)) +
# #   geom_line() + facet_wrap(~variable, scales='free') + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Mean concentration of fecal SCFAs over time')
# #  
# 
# 
# #######
# 
# D0vfas <- res.all %>% 
#   filter(time == 0 & treatment %in% c('control', 'RPS', 'Acid','RCS')) %>% 
#   mutate(treatment=factor(treatment, levels = c('control', 'RPS', 'Acid','RCS')))
# 
# # D0vfas$treatment
# 
# D21vfas <- vfas %>%
#   filter(hour == 0) %>% 
#   mutate(treatment=factor(treatment, levels = c('control', 'RPS', 'Acid','RCS')))
# 
# # D21vfas$treatment
# 
# D0_results <- D0vfas %>%
#   # select(-formate) %>%
#   filter(pignum != 101) %>% 
#   pivot_longer(cols = acetate:total,
#                names_to = 'VFA',
#                values_to = 'concentration') %>% 
#   group_by(VFA) %>% 
#   nest() %>% ungroup() %>% 
#   mutate(day=0, 
#          ANOVA=map(.x=data, ~ aov(data = .x, formula = concentration ~ treatment)),
#          TUK=map(.x = ANOVA, TukeyHSD), 
#          tid_TUK=map(.x=TUK, broom::tidy)) %>% 
#   select(-data, -ANOVA, -TUK) %>% 
#   unnest(cols = tid_TUK) %>% 
#   filter(grepl('control', contrast))
# 
# D21_results <- D21vfas %>% 
#   pivot_longer(cols = acetate:total,
#                names_to = 'VFA',
#                values_to = 'concentration') %>% 
#   group_by(VFA) %>% 
#   nest() %>% ungroup() %>% 
#   mutate(day=21, 
#          ANOVA=map(.x=data, ~ aov(data = .x, formula = concentration ~ treatment)), 
#          TUK=map(.x = ANOVA, TukeyHSD), 
#          tid_TUK=map(.x=TUK, broom::tidy)) %>% 
#   select(-data, -ANOVA, -TUK) %>% 
#   unnest(cols = tid_TUK) %>% 
#   filter(grepl('control', contrast))
# 
# 
# 
# D21_results %>% filter(VFA %in% c('butyrate', 'caproate', 'valerate')) %>% 
#   ggplot(aes(x=contrast, y=estimate, ymin=conf.low, ymax=conf.high, color=contrast)) + 
#   geom_hline(yintercept = 0, color='grey')+
#   geom_pointrange() + facet_wrap(~VFA, scales = 'free') + coord_flip()
# 
# 
# 
# # probably omit D0 results, maybe supplement if anything.  
# 
# 
# D0_results %>%# filter(VFA %in% c('butyrate', 'caproate', 'valerate')) %>% 
#   ggplot(aes(x=contrast, y=estimate, ymin=conf.low, ymax=conf.high, color=contrast)) + 
#   geom_hline(yintercept = 0, color='grey')+
#   geom_pointrange() + facet_wrap(~VFA, scales = 'free') + coord_flip()
# 
# 
# 
# 
# 
# D0_results %>% group_by(VFA) %>% 
#   summarise(MI=min(conf.low) , 
#             MA=max(conf.high), 
#             RANGE=abs(MI-MA), 
#             bump= RANGE*.1, 
#             MIN=MI - bump, 
#             MAX=MA + bump, 
#             LIM = ifelse(abs(MIN) > MAX, abs(MIN), abs(MAX)))
# 
