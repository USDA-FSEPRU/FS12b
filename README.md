# FS12b



esearch -db sra -query PRJNA638426| efetch -format runinfo > FS12b_SRA_meta.csv
mothur "#pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319, keepdots=F, processors=8);
                unique.seqs()"

