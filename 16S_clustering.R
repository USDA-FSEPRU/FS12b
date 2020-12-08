

clust_assign <- function(phy, DAY, TISSUE, filt =NULL, DIST_METH = 'jaccard'){
  
  
  MET <- as(sample_data(phy), "data.frame")
  
  
  phy <- prune_samples(FS12b@sam_data$day == DAY & FS12b@sam_data$tissue == TISSUE, phy)
  phy <- rarefy_even_depth(phy)
  
  test <- phyloseq::distance(phy, method = DIST_METH) %>% as.matrix()
  
  
  
  rownames(test)
  
  long_dist <- 
    test %>%
    as.data.frame() %>% 
    rownames_to_column(var='fromID') %>% 
    gather(-fromID, key=toID, value='distance') %>%
    mutate(weight=1-distance)
  
  
  if(!is.null(filt)){
    filt_dist <- 
      long_dist %>% 
      filter(distance < filt & distance != 0)
    
  } else {
    
    filt_dist <- 
      long_dist %>% 
      filter(distance != 0)
    
    
  }
  
  
  
  
  hist(filt_dist$weight)
  
  
  library(igraph)
  
  g <- igraph::graph_from_data_frame(d=filt_dist, directed = FALSE)
  
  
  # plot(g)
  
  # V(g)
  # 
  # E(g)
  # igraph::is.weighted(g)
  #
  clouv <- cluster_louvain(g)
  membership(clouv)
  
  clusts <- 
    tibble(ID=names(membership(clouv)), 
           clust=as.character(membership(clouv)))
  
  # MET$ID == names(membership(clouv))
  
  ###
  
  
  MET <- 
    MET %>% right_join(clusts)
  
  p1 <- 
    MET %>% group_by(day, clust) %>% 
    summarise(mean_AULC=mean(AULC), 
              min_AULC = min(AULC), 
              max_AULC = max(AULC)) %>% 
    ggplot(aes(x=clust, y=mean_AULC, ymin=min_AULC, ymax=max_AULC)) + 
    geom_pointrange(aes(color=clust)) 
  
  p2 <- 
    MET %>% group_by(day, clust, treatment) %>% 
    summarise(mean_AULC=mean(AULC), 
              min_AULC = min(AULC), 
              max_AULC = max(AULC), 
              stder    =std_err(AULC), 
              numobs   =n()) %>% 
    ggplot(aes(x=clust, y=mean_AULC, ymin=mean_AULC-stder, ymax=mean_AULC+stder, color=treatment)) + 
    geom_pointrange(aes(color=treatment)) + 
    geom_text(aes(label=numobs), nudge_y = 5)

    clust_frame <- MET %>% select(pignum, treatment, clust, AULC)
  
  return(list(p1,p2, clust_frame))
  
  
}




clust_assign(phy = FS12b, DAY = 'D0', TISSUE ='F')
clust_assign(phy = FS12b, DAY = 'D2', TISSUE ='F')
clust_assign(phy = FS12b, DAY = 'D7', TISSUE ='F')
clust_assign(phy = FS12b, DAY = 'D14', TISSUE ='F')
clust_assign(phy = FS12b, DAY = 'D21', TISSUE ='F')


clust_assign(phy = FS12b, DAY = 'D21', TISSUE ='X')
clust_assign(phy = FS12b, DAY = 'D21', TISSUE ='C')
clust_assign(phy = FS12b, DAY = 'D21', TISSUE ='I')


MET <- read_tsv('./data/FS12b_meta.tsv') %>%
  mutate(ID=sample_ID) %>%
  select(ID, everything()) %>% 
  mutate(treatment=factor(treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))) %>% 
  column_to_rownames(var='sample_ID') %>%
  filter(tissue =='F')



FECES_ONLY <- prune_samples(FS12b@sam_data$tissue == 'F', FS12b)


test <- phyloseq::distance(FECES_ONLY, method = 'jaccard') %>% as.matrix()

rownames(test)

long_dist <- 
  test %>%
  as.data.frame() %>% 
  rownames_to_column(var='fromID') %>% 
  gather(-fromID, key=toID, value='distance') %>%
  mutate(weight=1-distance)

filt_dist <- 
  long_dist #%>% 
  filter(distance < 0.5 & distance != 0)


hist(filt_dist$weight)


library(igraph)

g <- igraph::graph_from_data_frame(d=filt_dist, directed = FALSE)


plot(g)

V(g)

E(g)
igraph::is.weighted(g)
#
clouv <- cluster_louvain(g)
membership(clouv)

clusts <- 
  tibble(ID=names(membership(clouv)), 
         clust=as.character(membership(clouv)))

# MET$ID == names(membership(clouv))

###
library(Rtsne)

tsne <- Rtsne(X = FS12b@otu_table)

plot(tsne$Y)

MET <- 
  MET %>% left_join(clusts)


MET %>% group_by(day, clust) %>% filter(day != 'D0') %>% 
  summarise(mean_log_sal=mean(AULC)) %>% 
  ggplot(aes(x=clust, y=mean_log_sal)) + geom_point(aes(color=day)) 


MET %>% filter(clust == '4')
MET %>% filter(pignum == '97')








##############
#######EXAMPLE #########
# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(mydata, method.hclust="ward",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

