# FS12b


##### to fecth raw data SRA info
esearch -db sra -query PRJNA638426| efetch -format runinfo > FS12b_SRA_meta.csv  

  
##### command used to create silva reference alignment
mothur "#pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319, keepdots=F, processors=8);
                unique.seqs()"  
  

##### need to copy over data from mothur out
 1029  cp './raw_data/NEW_MOTHUR_OUT/FS12b.shared' ./data/
 1030  cp './raw_data/NEW_MOTHUR_OUT/FS12b_OTU.taxonomy' ./data/

