library(cowplot)
library(tidyverse)
pens <- read_csv('FS12_all_treatments.csv') %>% select(pen, treatment) %>% unique() %>% arrange(pen)

unique(pens$treatment)

pens <- 
  pens %>%
  mutate(treatment=ifelse(treatment %in% c('Control', 'RPS', 'Acid', 'RCS'), treatment, 'other'), 
         xpos=rep(c(1,2,3,4,5,6,7,8), each=12), 
         ypos=rep(c(1,2,3,4,5,6,7,8,9,10,11,12), times=8), 
         treatment = factor(treatment ,levels = c('Control', 'RPS', 'Acid', 'RCS', 'other')))


pens %>% mutate(cohort=ifelse(treatment != 'other', TRUE, FALSE))

pigpen <- 
  tibble(
  pignum=1:(96*5), 
  pen=rep(1:96, each=5), 
  cohort=rep(c(F,F,T,F,F), times=96), 
  pignumpen=rep(1:5, times=96),
  pos_adj=rep(c(-.3,-.15,0,.15,.3), times=96)) %>% 
  left_join(pens) %>% 
  mutate(xpigpos=xpos + pos_adj) %>% 
  mutate(cohort=ifelse(treatment == 'other', F, cohort), 
         strk=ifelse(cohort, 1.5,1))
pigpen

rep(c(1,1,1.5,1,1), times=96)

p1 <- 
  pens %>%
  ggplot() +
  geom_tile(aes(x=xpos, y=ypos, fill=treatment),width=.75, height=.75, color='black') + 
  geom_point(data=pigpen, aes(x=xpigpos, y=ypos, shape=cohort, stroke=strk))+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ylab('')+ xlab('') + theme_cowplot()+
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  annotate('text', x=4.5, y=15, label='Nursery Pen Layout', size=5)+
  annotate('text', x=4.5, y=13.75, label='Five pigs per pen, in nursery for 4 weeks', size=5)
  # ggtitle('Nursery Pen Layout', 'five pigs per pen, in Nursery for 4 weeks')

p1




#####

room_tibble <- 
  tibble(
    room=c(1:4), 
    xpos=c(0,2,4,6), 
    ypos=c(0,0,0,0), 
    treatment = factor(c('Control', 'RPS', 'Acid', 'RCS') ,levels = c('Control', 'RPS', 'Acid', 'RCS'))
  )

room_tibble

pig_room_tibble <- 
  tibble(
    room=rep(1:4, each=10), 
    xpos_adj=rep(c(-.5,-.5,-.25,-.25,0,0,.25,.25,.5,.5), times=4),
    ypos_adj=rep(c(rep(c(.25,-.25), times=20))),
    pignum_room=rep(1:10, times=4)
  ) %>% left_join(room_tibble) %>% 
  mutate(xpos=xpos+xpos_adj, 
         ypos=ypos+ypos_adj)

library(ggtext)


p2 <- 
  pig_room_tibble %>%
  ggplot()+
  geom_tile(data=room_tibble, aes(x=xpos, y=ypos, fill=treatment),color='black', width=1.5, height=.8, size=1)+
  geom_point(aes(x=xpos, y=ypos,fill=treatment), shape=24, size=3, stroke=1) +
  geom_text(data = room_tibble, aes(label=treatment, y=ypos, x=xpos), size=5)+
  annotate(geom='text', x=3, y=.6,size=5, label='One pig per pen from each treatment transfered to isolation rooms')+
  annotate(geom='richtext', x=3, y=.5,size=5, label='Challenged with 8e7 CFU <i>S. enterica</i> strain SX240', fill = NA, label.color = NA,)+
  annotate(geom='text', x=3, y=-.5,size=5, label='Feces collected at 0, 2, 7, 14, and 21 dpi')+
  annotate(geom='text', x=3, y=-.6,size=5, label='Necropsies performed at 21 dpi')+
  
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + 
  theme_cowplot() +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = 'none')
  # xlim(-1,7)+
  # ylim(-1,1) +
p2
# label="P(italic(i))==8~italic(i)", parse=TRUE



mfig <- ggdraw()+
  draw_plot(p1, 0,.4,1,.6)+
  draw_plot(p2, 0,0,1,.45)+
  draw_plot_label(x=c(0,0), y=c(1,.45), label = c('A', 'B'))
mfig





ggsave(mfig,
       filename = './output/figure1.jpeg',
       width = 180,
       height = 180,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

