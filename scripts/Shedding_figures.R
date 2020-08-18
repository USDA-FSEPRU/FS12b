require(pracma)
library(tidyverse)
library(Hmisc)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
# Read in data

sal_data <- read_csv('./data/FS12b_salmonella_data.csv')

# orders the treatment factors and makes a timepoint factor for plotting ease

sal_data$treatment <- factor(sal_data$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

sal_data$treatXtime <- factor(sal_data$treatXtime, levels = c('control_0', 'RPS_0', 'Acid_0', 'Zn+Cu_0', 'RCS_0', 'Bglu_0',
                                                              'control_2', 'RPS_2', 'Acid_2', 'Zn+Cu_2', 'RCS_2', 'Bglu_2',
                                                              'control_7', 'RPS_7', 'Acid_7', 'Zn+Cu_7', 'RCS_7', 'Bglu_7',
                                                              'control_14', 'RPS_14', 'Acid_14', 'Zn+Cu_14', 'RCS_14', 'Bglu_14',
                                                              'control_21', 'RPS_21', 'Acid_21', 'Zn+Cu_21', 'RCS_21', 'Bglu_21'))

sal_data$time_point_fact <- factor(sal_data$time_point)
# sal_data$time_point
#

sal_data <- sal_data %>% filter(!(treatment %in% c('Zn+Cu', 'Bglu')) & pignum != 101)




all_daily <- sal_data %>% filter(pignum != 101) %>% group_by(time_point, treatment) %>%
  summarise(mean_sal=mean(log_sal),
            sd_sal=sd(log_sal),
            num=n(),
            se_sal=sd_sal/sqrt(num))






### Figure 1A
F1A <- all_daily %>% 
  mutate(TPP=case_when(
    treatment == 'RPS'  ~    time_point, 
    treatment == 'Acid' ~    time_point - .15, 
    treatment == 'RCS'  ~    time_point + .15, 
    treatment == 'control' ~ time_point, 
    TRUE ~ 1000000000000000000)) %>% 
  filter(!(treatment %in% c('Zn+Cu', 'Bglu'))) %>%
  ggplot(aes(x=TPP, y=mean_sal, fill=treatment, group=treatment)) +
  # geom_jitter(aes(x=time_point, y=log_sal, color=treatment), data=filter(sal_data, !(treatment %in% c('Zn+Cu', 'Bglu'))& pignum !=101), alpha=.5) + 
  geom_line(aes(color=treatment), size=1.5)+
  geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) +
  geom_point(shape=21, size=4) +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ylab('log CFU / g feces') +
  xlab('Day post-challenge') +
  theme(legend.position = 'none', 
        axis.title = element_text(size = 11), 
        panel.grid.major = element_line(color='grey'))
  # ggtitle('Daily shedding, group summary statistics')


F1A
###


#############






#####
# use this #
# nest and do an anova at each timepoint
library(broom)
daily_tests <- sal_data %>% filter(time_point != 0) %>% group_by(time_point) %>% nest() %>% 
  mutate(AOV=map(.x = data, ~ aov(data=.x, log_sal ~ treatment)), 
         tid_AOV=map(AOV, tidy), 
         TUK=map(AOV, TukeyHSD), 
         tid_TUK=map(TUK, tidy))


daily_tests %>% select(time_point, tid_AOV) %>% unnest(cols = tid_AOV)



daily_tuks <- daily_tests %>%
  select(time_point, tid_TUK) %>% unnest(cols = tid_TUK) %>%
  filter(grepl('control', comparison)) %>% 
  mutate(tuk_pval=adj.p.value, 
         fdr_pval=p.adjust(adj.p.value, method = 'fdr'), 
         TPP=paste('Day', time_point), 
         TPP=factor(TPP, levels = c('Day 2', 'Day 7', 'Day 14', 'Day 21'))) %>% 
  select(-adj.p.value)

anno <- tibble(y=c(0,0,0,0), 
               TPP=factor(c('Day 2','Day 7','Day 14','Day 21'), levels = c('Day 2', 'Day 7', 'Day 14', 'Day 21')), 
               x=c(3.6,3.6,3.6,3.6), 
               labtext = c('ANOVA p=0.10', 
                           'ANOVA p=0.07',
                           'ANOVA p=0.26',
                           'ANOVA p=0.11'))

### FIG 1B
F1B <- daily_tuks %>%
  ggplot(aes(x=comparison, y=estimate, ymin=conf.low, ymax=conf.high, color=comparison)) +
  geom_hline(yintercept = 0, color='grey')+
  geom_pointrange(size=.5) + 
  geom_label(data=anno, aes(x=x,y=y,label=labtext), inherit.aes = FALSE)+
  geom_text(aes(label=round(tuk_pval, digits = 2), y=2))+
  coord_flip() + scale_x_discrete(expand = expand_scale(add = c(.5,1)))+
  facet_wrap(.~TPP, ncol = 4) + 
  ylim(-3,3) + scale_color_manual(values=c('red','orange','#3399FF')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = 'none', 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size=11), 
        panel.grid.major.x = element_line(color = 'grey', size = .2))


F1B


# FIGURE 2A


sum_sal <- sal_data %>% group_by(pignum) %>%
  summarise(AULC=trapz(time_point, log_sal),
            sum=sum(log_sal), 
            maxshed=log_sal[which.max(log_sal)], 
            day_max=time_point[which.max(log_sal)], 
            pos_samples=sum(log_sal > 0), 
            treatment=unique(treatment))


# re orders sum_sal
sum_sal <- sum_sal[match(filter(sal_data, time_point==2)$pignum,sum_sal$pignum),]



sum_sal$treatment <- factor(sum_sal$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))


F2A <- filter(sum_sal, pignum !=101) %>% 
  ggplot(aes(x=treatment, y=AULC, fill=treatment))+
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(fill=treatment), shape=21, size=1.75, stroke=1, width = .13) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  # ggtitle('Cumulative Salmonella shedding (AULC)', subtitle = 'ANOVA P = 0.0123') +
  theme(legend.position = 'none', 
        panel.border = element_rect(colour = 'black', fill=NA), 
        axis.title.x = element_blank(), 
        panel.grid.major = element_line(color='grey'))


F2A


### USE THIS FOR AULC STAT TEST ###
aov_AULC <- aov(data=sum_sal, AULC~treatment)
summary(aov_AULC)

AULC_tuk <- TukeyHSD(aov_AULC) %>% tidy()

F2B <- AULC_tuk %>% filter(grepl('control', comparison)) %>% 
  mutate(comparison=factor(comparison, levels = c('RCS-control', 'Acid-control','RPS-control'))) %>%
  ggplot(aes(x=comparison, y=estimate, ymin=conf.low, ymax=conf.high,color=comparison)) +
  geom_hline(yintercept = 0, color='grey') +
  geom_pointrange(size=.75) +
  geom_text(aes(x=comparison,y=8,
                label=paste('P=', round(adj.p.value, 2))), 
            hjust=0) +
  ggtitle('ANOVA P = 0.012') + 
  ylim(-32,32)+
  # scale_y_continuous(expand = expand_scale(add = c(5,15)))+
  theme(panel.border = element_rect(color='black', fill=NA), 
        legend.position='none', 
        axis.title.y = element_blank(), 
        axis.text = element_text(size=10), 
        panel.grid.major.x = element_line(color='grey', size=.2))+
  coord_flip() +
  scale_color_manual(values=c('red','orange','#3399FF')) 
F2B







################################
# setwd("~/FS12")
# TISSUES


tis <- read.csv('./data/D21_tissues.csv')
tis$log_sal <- log10(tis$Salmonella)
tis$log_sal <- as.numeric(sub(-Inf, 0, tis$log_sal))
tis <- na.exclude(tis)

tis <- tis %>% filter(!treatment %in% c('Zn+Cu', 'Bglu'))

tis$treatment <- factor(tis$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))


# pairwise.wilcox.test(tis$log_sal[tis$tissue == 'cecal_cont'], tis$treatment[tis$tissue == 'cecal_cont'], p.adjust.method = 'none')
# 
# pairwise.t.test(tis$log_sal[tis$tissue == 'cecal_cont'], tis$treatment[tis$tissue == 'cecal_cont'], p.adjust.method = 'none')
# pairwise.t.test(tis$log_sal[tis$tissue == 'Cecum'], tis$treatment[tis$tissue == 'Cecum'], p.adjust.method = 'none')




# tis %>% filter(tissue=='cecal_cont') %>%
#   ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Cecal Contents')
# 

# tis %>% #filter(tissue=='cecal_cont') %>%
#   ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +
#   geom_jitter(shape=21,width = .2, size=2.25) +
#   facet_wrap(~tissue) +
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Salmonella colonization at D21, tissues')
# 

# tables?

tissum <- tis %>% group_by(treatment, tissue) %>% 
  summarise(mean_lcfu=mean(log_sal),
            num_pos=sum(log_sal>0), 
            tot_obs=n())


tissum %>% filter(tissue=='cecal_cont')




#

tis_tests <- tis %>%  group_by(tissue) %>% nest() %>% 
  mutate(AOV=map(.x = data, ~ aov(data=.x, log_sal ~ treatment)), 
         tid_AOV=map(AOV, tidy), 
         TUK=map(AOV, TukeyHSD), 
         tid_TUK=map(TUK, tidy))


tis_tests %>% select(tissue, tid_AOV) %>% unnest(cols = tid_AOV)

tis_anno <- tibble(y=c(0,0,0,0,0), 
                   tissue=factor(c('cecal_cont','Cecum','ICLN','IPP', 'Tonsil'), levels = c('cecal_cont','Cecum','ICLN','IPP', 'Tonsil')), 
                   x=c(3.9,3.9,3.9,3.9,3.9), 
                   labtext = c('ANOVA p=0.02', 
                               'ANOVA p=0.05',
                               'ANOVA p=0.13',
                               'ANOVA p=0.54',
                               'ANOVA p=0.01'))


  
tissue_tuks <- tis_tests %>%
  select(tissue, tid_TUK) %>% unnest(cols = tid_TUK) %>%
  filter(grepl('control', comparison)) %>% 
  mutate(tuk_pval=adj.p.value, 
         fdr_pval=p.adjust(adj.p.value, method = 'fdr')) %>% 
  select(-adj.p.value)

tissue_tuks
  
#

all_tis <- tis %>% group_by(tissue, treatment) %>%
  summarise(mean_sal=mean(log_sal),
            sd_sal=sd(log_sal),
            num=n(),
            se_sal=sd_sal/sqrt(num))






### Figure 3A
F3A <- all_tis %>%
  ggplot(aes(x=treatment, y=mean_sal, fill=treatment, group=treatment)) +
  # geom_jitter(aes(x=time_point, y=log_sal, color=treatment), data=filter(sal_data, !(treatment %in% c('Zn+Cu', 'Bglu'))& pignum !=101), alpha=.5) + 
  geom_col(aes(fill=treatment), size=1.5)+
  geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ylab('log CFU / g tissue') +
  theme(panel.border = element_rect(color = 'black', fill = NA), 
        axis.text.x = element_text(angle = -45, hjust = 0), 
        axis.title.x = element_blank(), 
        legend.position = 'none', 
        axis.title.y = element_text(size = 11), 
        panel.grid.major = element_line(color='grey', size=.25))+
  # ggtitle('tissue colonization') +
  facet_wrap(~tissue, ncol = 5)


F3A



#
# Figure 3B

F3B <- tissue_tuks %>%
  mutate(comparison = factor(comparison, levels=c('RCS-control', 'Acid-control','RPS-control'))) %>% 
  ggplot(aes(x=comparison, y=estimate, ymin=conf.low, ymax=conf.high, color=comparison)) +
  geom_hline(yintercept = 0, color='grey')+
  geom_pointrange(size=.5) + 
  geom_label(data=tis_anno, aes(x=x,y=y,label=labtext), inherit.aes = FALSE, size=3)+
  geom_text(aes(label=round(tuk_pval, digits = 2), y=2))+
  coord_flip() + scale_x_discrete(expand = expand_scale(add = c(1,1.5)))+
  facet_wrap(.~tissue, ncol = 5) +
  ylim(-3.5,3.5) + scale_color_manual(values=c('red','orange','#3399FF')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y=element_blank(), 
        legend.position = 'none', 
        panel.grid.major = element_line(color='grey', size=.25))


F3B

### COWPLOT ZONE ###

# Figure 1

# F1A
# F1B


fig_1 <- ggdraw()+
  draw_plot(F1A, 0,.45,1,.55)+
  draw_plot(F1B, 0,0,1,.45)+
  draw_plot_label(x=c(0,0), y=c(1,.45), label = c('A', 'B'))
fig_1


ggsave(fig_1,
       filename = './output/figs/figure1.jpeg',
       width = 180,
       height = 150,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')


# Figure 2

# F2A
# F2B


fig_2 <- ggdraw()+
  draw_plot(F2A, 0,0,.6,1)+
  draw_plot(F2B, .6,0,.4,1)+
  draw_plot_label(x=c(0,.6), y=c(1,1), label = c('A', 'B'))
fig_2


ggsave(fig_2,
       filename = './output/figs/figure2.jpeg',
       width = 180,
       height = 85,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')



# Figure 3

# F3A
# F3B


fig_3 <- ggdraw()+
  draw_plot(F3A, 0,.5,1,.5)+
  draw_plot(F3B, 0,0,1,.5)+
  draw_plot_label(x=c(0,0), y=c(1,.5), label = c('A', 'B'))
fig_3


ggsave(fig_3,
       filename = './output/figs/figure3.jpeg',
       width = 180,
       height = 150,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')



########### FOR A HIGH VS LOW SHEDDER SHED FIGURE ##########

tis_RPS <- tis %>% filter(treatment %in% c('control', 'RPS'))
sum_sal_RPS <- sum_sal %>% filter(treatment %in% c('control', 'RPS'))

controls <- sum_sal$pignum[sum_sal_RPS$treatment == 'control']


tis_RPS$shed <- ifelse(tis_RPS$pignum %in% c(373,321,181,392,97), 'RPS_low',
                       ifelse(tis_RPS$pignum %in% controls,'control', 'RPS_high'))
                       
                       
sum_sal_RPS$shed <- ifelse(sum_sal_RPS$pignum %in% c(373,321,181,392,97), 'RPS_low',
                           ifelse(sum_sal_RPS$pignum %in% controls,'control', 'RPS_high'))

tis_RPS %>% ggplot(aes(x=shed, y=log_sal, group=shed, fill=shed)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .1, size=2.25) +
  facet_wrap(~tissue) +
  scale_fill_manual(values=c('#33CC33', '#246DB6', '#47D6FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Cecal Contents')

sum_sal_RPS %>% ggplot(aes(x=shed, y=AULC, group=shed, fill=shed)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .1, size=2.25) +
  #facet_wrap(~tissue) +
  scale_fill_manual(values=c('#33CC33', '#246DB6', '#47D6FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('AULC')


#C8C226
#985CF9

# CONTROL COLOR    '#33CC33'
# RPS COLOR        '#3399FF'
# HIGH SHED COLOR: '#246DB6'
# LOW SHED COLOR:  '#47D6FF'

darken('#3399FF')
lighten('#3399FF')



sum_sal %>% ggplot(aes(x=0, y=AULC,fill=treatment)) +
  # geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .01, size=2.25)+ geom_violin()



sum_sal %>% 
  ggplot(aes(x=AULC, fill=treatment)) +
  geom_histogram(fill='grey', color='grey') +
  geom_histogram(color=alpha('black', alpha = .5)) +
  # geom_vline(xintercept = c(30.75, 44.5), color='purple') +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))



# 
# 
# ### NEED CECAL VFA DATA FOR THIS NEXT SECTION ###
# 
# vfas_for_cor <- read_csv('data/FS12b_vfas_for_cor.csv')
# 
# ### NEED TO REMOVE ONE TREATMENT BEFORE MERGE
# vfas_for_cor <- vfas_for_cor %>% select(-treatment)
# sal_for_cor <- merge(sal_for_cor, vfas_for_cor, by = 'pignum')
# rownames(sal_for_cor) <- sal_for_cor$pignum
# 
# 
# ##### fecal corrs #####
# fec_cor <- res.all %>% select(-treatment) %>% full_join(sum_sal, by = 'pignum')
# fec_cor <- fec_cor[-288,]
# 
# sum_sal
# 
# fec_cor %>% filter(time == 0) %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
# fec_cor %>% filter(time == 0) %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
# fec_cor %>% filter(time == 0) %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
# 
# 
# fec_cor %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + facet_wrap(~time)
# fec_cor %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + facet_wrap(~time)
# 
# 
# 
# fec_cor %>% filter(time == 0) %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
# 
# 
# ### is this meta from the 16S stuff?
# meta %>% filter(experiment == 'X12b') %>% ggplot(aes(x=caproate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day)
# 
# meta_for_corr <- meta %>% mutate(day_fact=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% filter(experiment == 'X12b')
# 
# meta_for_corr %>% filter(tissue =='F' & day != 'D0') %>% ggplot(aes(x=caproate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day_fact)
# meta_for_corr %>% filter(tissue =='F' & day != 'D0') %>% ggplot(aes(x=valerate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day_fact)
# meta_for_corr %>% filter(tissue =='F' & day != 'D0') %>% ggplot(aes(x=butyrate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day_fact)
# 
# # meta_for_corr$treatment
# sal_vfa_cor <- meta_for_corr %>% select(pignum,day, tissue, treatment, log_sal, AULC, ends_with('ate')) %>%
#   filter(pignum!=101 & tissue != 'Q' & treatment %in% c('Control', 'RPS', 'Acid', 'RCS')) %>% na.omit()
# 
# sal_vfa_cor <- sal_vfa_cor %>% mutate(total=rowSums(.[grep("ate", names(.))]))
# 
# 
# globa_log_sal_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue) %>% 
#                                           summarise(ace_corP=cor.test(log_sal, acetate)$p.value, 
#                                                     pro_corP=cor.test(log_sal, propionate)$p.value, 
#                                                     but_corP=cor.test(log_sal, butyrate)$p.value, 
#                                                     val_corP=cor.test(log_sal, valerate)$p.value, 
#                                                     cap_corP=cor.test(log_sal, caproate)$p.value, 
#                                                     isob_corP=cor.test(log_sal, isobutyrate)$p.value, 
#                                                     isov_corP=cor.test(log_sal, isovalerate)$p.value, 
#                                                     tot_corP=cor.test(log_sal, total)$p.value)%>% na.omit() %>%
#   gather(-(day:tissue), key = 'vfa', value = 'pval') %>% filter(pval < 0.05)
# 
# 
# 
# globa_AULC_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue) %>% 
#   summarise(ace_corP=cor.test(AULC, acetate, method = 'spearman')$p.value, 
#             pro_corP=cor.test(AULC, propionate, method = 'spearman')$p.value, 
#             but_corP=cor.test(AULC, butyrate, method = 'spearman')$p.value, 
#             val_corP=cor.test(AULC, valerate, method = 'spearman')$p.value, 
#             cap_corP=cor.test(AULC, caproate, method = 'spearman')$p.value, 
#             isob_corP=cor.test(AULC, isobutyrate, method = 'spearman')$p.value, 
#             isov_corP=cor.test(AULC, isovalerate, method = 'spearman')$p.value, 
#             tot_corP=cor.test(AULC, total)$p.value, method = 'spearman')%>% na.omit() %>%
#   gather(-(day:tissue), key = 'vfa', value = 'pval') %>% filter(pval < 0.05)
# 
#  c('acetate', 'propionate', 'valerate', 'caproate', 'total')
#  
# globa_AULC_VFA_corrs %>% filter(day =='D0')
# globa_AULC_VFA_corrs %>% filter(day =='D2')
# 
# 
# treat_log_sal_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue, treatment) %>% 
#                                           summarise(ace_corP=cor.test(log_sal, acetate, method = 'spearman')$p.value, 
#                                                     pro_corP=cor.test(log_sal, propionate, method = 'spearman')$p.value, 
#                                                     but_corP=cor.test(log_sal, butyrate, method = 'spearman')$p.value, 
#                                                     val_corP=cor.test(log_sal, valerate, method = 'spearman')$p.value, 
#                                                     cap_corP=cor.test(log_sal, caproate, method = 'spearman')$p.value, 
#                                                     isob_corP=cor.test(log_sal, isobutyrate, method = 'spearman')$p.value, 
#                                                     isov_corP=cor.test(log_sal, isovalerate, method = 'spearman')$p.value, 
#                                                     tot_corP=cor.test(log_sal, total, method = 'spearman')$p.value) %>% na.omit() %>%
#   gather(-(day:treatment), key = 'vfa', value = 'pval') %>% filter(pval < 0.15)
# 
# 
# treat_AULC_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue, treatment) %>% 
#   summarise(ace_corP=cor.test(AULC, acetate, method = 'spearman')$p.value, 
#             pro_corP=cor.test(AULC, propionate, method = 'spearman')$p.value, 
#             but_corP=cor.test(AULC, butyrate, method = 'spearman')$p.value, 
#             val_corP=cor.test(AULC, valerate, method = 'spearman')$p.value, 
#             cap_corP=cor.test(AULC, caproate, method = 'spearman')$p.value, 
#             isob_corP=cor.test(AULC, isobutyrate, method = 'spearman')$p.value, 
#             isov_corP=cor.test(AULC, isovalerate, method = 'spearman')$p.value, 
#             tot_corP=cor.test(AULC, total, method = 'spearman')$p.value) %>% na.omit() %>%
#   gather(-(day:treatment), key = 'vfa', value = 'pval') %>% filter(pval < 0.15)
# 
# ##### schemezone #####
# globa_AULC_VFA_corrs %>% filter(pval < 0.05 & day == 'D0')
# treat_AULC_VFA_corrs %>% filter(pval < 0.15 & day == 'D0')
# 
# globa_AULC_VFA_corrs %>% filter(pval < 0.05 & day=='D21') %>% print(n=40)
# treat_AULC_VFA_corrs %>% filter(pval < 0.15 & day=='D21') %>% print(n=40)
# 
# treat_log_sal_VFA_corrs %>%filter(pval <0.1) %>%  print(n = 40)
# globa_log_sal_VFA_corrs
# 
# #### AULC CORRELATIONS PLOTS #####
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D0') %>% 
#   ggplot(aes(x=mM, y=AULC)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +geom_point(aes(color=treatment), alpha=.5)+
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D0 fecal VFAs correlate with final AULC')
# 
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='C') %>% 
#   ggplot(aes(x=mM, y=AULC)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +geom_point(aes(color=treatment), alpha=.5)+
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D21 cecal VFAs correlate with final AULC')
# 
# #### LOG SAL CORRELATIONS PLOTS #########
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D2') %>% 
#   ggplot(aes(x=mM, y=log_sal)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D2 fecal VFAs correlate with log_sal')
# 
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D7') %>% 
#   ggplot(aes(x=mM, y=log_sal)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D7 fecal VFAs correlate with log_sal')
# 
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D14') %>% 
#   ggplot(aes(x=mM, y=log_sal)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D14 fecal VFAs correlate with log_sal')
# 
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='F') %>% 
#   ggplot(aes(x=mM, y=log_sal)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D21 fecal VFAs correlate with log_sal')
# 
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='X') %>% 
#   ggplot(aes(x=mM, y=log_sal)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D21 cecal VFAs correlate with cecal_tissue log_sal')
# 
# sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='C') %>% 
#   ggplot(aes(x=mM, y=log_sal)) +
#   geom_smooth(color='black', method = 'lm', fill=NA) +
#   geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
#   facet_wrap(~VFA, scales = 'free') + ggtitle('D21 cecal VFAs correlate with cecal log_sal')
# 
# 
# 
# ########
# treat_log_sal_VFA_corrs %>% na.omit()
# 
# ggplot(sal_for_cor, aes(x=butyrate, y=AULC)) +
#   geom_point(aes(color=treatment)) + geom_smooth(method = 'lm') + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Correlation between cecal butyrate concentration and total Salmonella shedding', subtitle = 'Spearman: -0.4, p=0.003')
# 
# 
# ggplot(sal_for_cor, aes(x=caproate, y=AULC)) +
#   geom_point(aes(color=treatment)) + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Correlation between cecal caproate concentration and total Salmonella shedding', subtitle = 'Spearman: -0.53, p=0.0005')
# ggplot(sal_for_cor, aes(x=valerate, y=AULC)) +
#   geom_point(aes(color=treatment)) + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Correlation between cecal valerate concentration and total Salmonella shedding', subtitle = 'Spearman: -0.38, p=0.0005')
# ggplot(sal_for_cor, aes(x=total, y=AULC)) +
#   geom_point(aes(color=treatment)) + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Correlation between cecal total SCFA concentration and total Salmonella shedding', subtitle = 'Spearman: -0.4, p=0.002')
# 
# ##### Shawn FSIS #######
# 
# ggplot(sal_for_cor, aes(x=butyrate, y=AULC)) +
#   geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme_bw(base_size = 16) + xlab('butyrate (mM)') + ggtitle('AULC vs Cecal butyrate (D21)')
# 
# ggplot(sal_for_cor, aes(x=total, y=AULC)) +
#   geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme_bw(base_size = 16) + xlab('Total SCFAs (mM)') + ggtitle('')
# 
# 
# 
# ggplot(sal_for_cor, aes(x=butyrate, y=AULC, color=treatment)) +
#   geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#     theme(axis.text = element_text(size=16), 
#         axis.title = element_text(size=16), 
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + ggtitle('') + xlab('butyrate (mM)')
# 
# ggplot(sal_for_cor, aes(x=valerate, y=AULC, color=treatment)) +
#   geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme(axis.text = element_text(size=16), 
#         axis.title = element_text(size=16), 
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + ggtitle('') + xlab('valerate (mM)')
# 
# ggplot(sal_for_cor, aes(x=caproate, y=AULC, color=treatment)) +
#   geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme(axis.text = element_text(size=16), 
#         axis.title = element_text(size=16), 
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + ggtitle('') + xlab('caproate (mM)')
# 
# ggplot(sal_for_cor, aes(x=total, y=AULC, color=treatment)) +
#   geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme(axis.text = element_text(size=16), 
#         axis.title = element_text(size=16), 
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + ggtitle('') + xlab('total (mM)')
# 
# 
# ############ TESTESTESTEST #################
# 
# 
# 
# sal_for_cor %>% filter(treatment %in% c('control','RPS')) %>% ggplot(aes(x=butyrate, y=AULC)) +
#   geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme_bw()+ ggtitle('Correlation between AULC and cecal butyrate at D21') + xlab('butyrate (mM)') 
# 
# 
# sal_for_cor %>% filter(treatment == 'Bglu') %>% ggplot(aes(x=butyrate, y=AULC)) +
#   geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme(axis.text = element_text(size=16), 
#         axis.title = element_text(size=16), 
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + ggtitle('') + xlab('butyrate (mM)')
# 
# 
# sal_for_cor %>% group_by(treatment) %>% summarise(butP = cor.test(AULC, butyrate)$p.value, 
#                                                   butT = cor.test(AULC, butyrate)$statistic,
#                                                   propP = cor.test(x=AULC, y=propionate)$p.value,
#                                                   propT = cor.test(x=AULC, y=propionate)$statistic,
#                                                   capP = cor.test(x=AULC, y=caproate)$p.value,
#                                                   capT = cor.test(x=AULC, y=caproate)$statistic,
#                                                   valP = cor.test(x=AULC, y=valerate)$p.value,
#                                                   valT = cor.test(x=AULC, y=valerate)$statistic,
#                                                   totP = cor.test(x=AULC, y=total)$p.value,
#                                                   totT = cor.test(x=AULC, y=total)$statistic,
#                                                   Pbutp = cor.test(x=AULC, y=P_butyrate)$p.value,
#                                                   PbutT = cor.test(x=AULC, y=P_butyrate)$statistic)
# 
# 
# 
# ggplot(sal_for_cor, aes(x=total, y=AULC)) +
#   geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   theme(axis.text = element_text(size=16), 
#         axis.title = element_text(size=16), 
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + ggtitle('') + xlab('Total SCFAs (mM)')
# 
# # filter(sal_for_cor, treatment %in% c('control', 'RPS'))
# # filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# # filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# 
# 
# # filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# # filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# # 
# 
# 
# # filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# # filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# # 
# 
# # filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# # filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# 
# 
# cor.test(filter(sal_for_cor, treatment =='control')$AULC,
#          filter(sal_for_cor, treatment =='control')$butyrate)
# 
# cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
#          filter(sal_for_cor, treatment =='RPS')$butyrate)
# 
# 
# cor.test(filter(sal_for_cor, treatment =='control')$AULC,
#          filter(sal_for_cor, treatment =='control')$valerate)
# 
# cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
#          filter(sal_for_cor, treatment =='RPS')$valerate)
# 
# cor.test(filter(sal_for_cor, treatment =='control')$AULC,
#          filter(sal_for_cor, treatment =='control')$caproate)
# 
# cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
#          filter(sal_for_cor, treatment =='RPS')$caproate)
# 
# cor.test(filter(sal_for_cor, treatment =='control')$AULC,
#          filter(sal_for_cor, treatment =='control')$total)
# 
# cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
#          filter(sal_for_cor, treatment =='RPS')$total)
# 
# 
# 
# cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
#          filter(sal_for_cor, treatment%in%c('RPS', 'control'))$total)
# 
# 
# 
# cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
#          filter(sal_for_cor, treatment%in%c('RPS', 'control'))$butyrate)
# 
# 
# 
# cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
#          filter(sal_for_cor, treatment%in%c('RPS', 'control'))$valerate)
# 
# 
# 
# cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
#          filter(sal_for_cor, treatment%in%c('RPS', 'control'))$caproate)
# 
# 
