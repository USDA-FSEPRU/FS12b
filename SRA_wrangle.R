library(tidyverse)


SRA_meta <-
  read_csv('./FS12b_SRA_meta.csv') %>% select(SampleName, Run)

FS12b_meta <- 
  read_tsv('./data/FS12b_meta.tsv') %>% 
  mutate(SampleName=sample_ID) %>% 
  left_join(SRA_meta) %>%
  write_tsv('./data/FS12b_meta.tsv')


FS12b_meta$Run %>% 
  write_lines('SRA_accessions.txt')

which(colnames(FS12b_meta)=='Run')


tibble(sample_ID=FS12b_meta$sample_ID,
       R1=paste0(FS12_meta$Run, '_1.fastq'),
       R2=paste0(FS12_meta$Run, '_2.fastq')) %>% 
  write_tsv('./raw_data/FS12b.files', col_names = FALSE)




