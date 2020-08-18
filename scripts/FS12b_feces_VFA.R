getwd()
# setwd('./FS12/FS12b/')

library(tidyverse)


####### WARNING! I THINK THERE IS A PROBLEM WITH THE DATA FROM DAY 7 and 21!!!!!

### SUSPICIOUS AMOUNTS OF LACTATE AND SUCCINATE ARE PRESENT, USUALLY THESE ARE VERY RARE
### I THINK THERE WAS AN ISSUE DURING SAMPLE PREP
### MY GUESS IS THAT THE SAMPLES WERE ALLOWED TO THAW AT ROOM TEMP FOR TOO LONG AND 
### SOME BACTERIAL METABOLIC ACTIVITY GENERATED EXCESS LACTATE AND SUCCINATE
### USUALLY THESE WOULD BE CONVERTED TO OTHER SCFAS BY ANAEROBIC FERMENTATION BUT
### THE ANAEROBES PROBABLY DIED DUE TO OXYGEN EXPOSURE.


# write_csv(res.all, 'FS12b_fecal_vfas.csv')

res.all <- read_csv('./data/FS12b_fecal_vfas.csv')
tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$butyrate, .$treatment, p.adjust.method = 'none'))
tests[[1]]


tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$caproate, .$treatment, p.adjust.method = 'none'))
#pairwise.wilcox.test(p.adjust.method = )
tests[[1]]


tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$valerate, .$treatment, p.adjust.method = 'none'))
#pairwise.wilcox.test(p.adjust.method = )
tests[[1]]


res.all$treatment <- factor(res.all$treatment, levels=c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))
res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=butyrate, fill=treatment)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + geom_jitter(shape=21, stroke=1.2, width = .2, size=2) + ggtitle("D0 Fecal Butyrate")
res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=caproate, fill=treatment)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+ geom_jitter(shape=21, stroke=1.2, width = .2, size=2) + ggtitle("D0 Fecal Caproate")
res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=valerate, fill=treatment)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+ geom_jitter(shape=21, stroke=1.2, width = .2, size=2) + ggtitle("D0 Fecal Valerate")


res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))


res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))

res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))

res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()


res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()

res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()


res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot() + geom_text(aes(label=pignum))
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))




res.all[res.all$treatment == 'asdfsa',]

res.all %>% filter(time == 2)
126 %in% res.all[res.all$time == 0, ]$pignum
126 %in% res.all[res.all$time == 2, ]$pignum
126 %in% res.all[res.all$time == 7, ]$pignum
126 %in% res.all[res.all$time == 14, ]$pignum
126 %in% res.all[res.all$time == 21, ]$pignum


326 %in% res.all[res.all$time == 0, ]$pignum
326 %in% res.all[res.all$time == 2, ]$pignum
326 %in% res.all[res.all$time == 7, ]$pignum
326 %in% res.all[res.all$time == 14, ]$pignum
326 %in% res.all[res.all$time == 21, ]$pignum

219 %in% res.all[res.all$time == 0, ]$pignum
219 %in% res.all[res.all$time == 2, ]$pignum
219 %in% res.all[res.all$time == 7, ]$pignum
219 %in% res.all[res.all$time == 14, ]$pignum
219 %in% res.all[res.all$time == 21, ]$pignum


392 %in% res.all[res.all$time == 0, ]$pignum
392 %in% res.all[res.all$time == 2, ]$pignum
392 %in% res.all[res.all$time == 7, ]$pignum
392 %in% res.all[res.all$time == 14, ]$pignum
392 %in% res.all[res.all$time == 21, ]$pignum

library(vegan)

res.all %>% ggplot(aes(x=treatment, y=acetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot() + geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=propionate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=butyrate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=valerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=caproate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=isovalerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=oxalate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=phenylacetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=succinate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=fumarate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=lactate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))


res.gather <- res.all %>% gather(key = variable, value = concentration, -(pignum:set))

res.gather %>% filter(variable == 'total')%>% ggplot(aes(x=time, y=concentration, group=pignum, color=treatment)) + geom_line()

res.sum <- res.gather %>%  group_by(time, treatment, variable) %>%
  summarise(mean=mean(concentration), sd=sd(concentration), n=n(), sterr=sd/sqrt(n))

res.sum <- res.sum %>% mutate(treat_var=paste(treatment, variable, sep='_'))

res.sum$treatment <- factor(res.sum$treatment, levels= c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

res.sum %>% ggplot(aes(x=time, y=mean, group=treat_var, color=treatment)) +
  geom_line() + facet_wrap(~variable, scales='free') + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Mean concentration of fecal SCFAs over time')
# 
# SCFA.j <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'jaccard')
# SCFA.j <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'time', distance_method = 'jaccard')
# 
# SCFA.j[[1]]$time <- factor(SCFA.j[[1]]$time)
# 
# 
# ggplot(SCFA.j[[1]], aes(x=MDS1, y=MDS2, color=treatment)) + geom_point()
# ggplot(SCFA.j[[1]], aes(x=MDS1, y=MDS2, color=time)) + geom_point() + geom_text(aes(label=pignum))
# 
# 
# SCFA.j[[2]]$treatment <- gsub('([0-9]+)_([A-Za-z]+)','\\2',SCFA.j[[2]]$group)
# SCFA.j[[2]]$time <- gsub('([0-9]+)_(.*)','\\1',SCFA.j[[2]]$group)
# ggplot(SCFA.j[[1]], aes(x=MDS1, y=MDS2, color=time)) + geom_point() +
#   geom_path(data = SCFA.j[[2]], aes(x=NMDS1, y=NMDS2, group=group))
# 
# 
# 
# ifelse(FS12b_HL@sam_data$pignum %in% c(373,321,181,392,97), 'low', 
#        ifelse(FS12b_HL@sam_data$pignum %in% c(50, 93,355, 244), 'high', 'Control'))


res.all$day <- paste('D',res.all$time, sep = '')
res.all$pig_day_tissue <- paste(res.all$pignum, res.all$day, 'F', sep = '_')
mergeme <- res.all %>% select(pig_day_tissue, ends_with('ate'))

