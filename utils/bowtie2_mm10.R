###:: mapping_mm10.sciTIP.R ::#####
####################################
### run after deindexing script
####################################

library("travis")
library("glue")
set.cores=20

# add folder containing FASTQ files
fastqs = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/fastq/deindexed_fastq/row_H"

# bowtie2 index
btindex.path="/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10_bowtie_index/mm10"
chromsize.path="/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt"

# system("mkdir -p ../data/mapped_mm10")
# system("mkdir -p ../data/stats_mapping")

fq = unique(strsplit(list.files(path = fastqs, pattern = "R1_001.fastq.gz", full.names = TRUE), '_R1_001.fastq.gz'))
length(fq) # number of cells

for(i in 1:length(fq)){
  Sam=bowtie2(read1files = paste0(fq[i],"_R1_001.fastq.gz"), read2files = paste0(fq[i],"_R2_001.fastq.gz"), indexfile = btindex.path, 
              alignMode = "very-sensitive-local", appendIndexToName=TRUE, unaligned =FALSE, 
              reorder=TRUE, threads=set.cores, minInsertSize = 10, maxInsertSize = 700)
  #system(glue("bowtie2 --threads 8 --very-sensitive-local --minins 10 --maxins 700 -x {btindex.path} -1 {fq[i]}_R1_001.fastq.gz -2 {fq[i]}_R2_001.fastq.gz -S {fq[i]}_mm10.sam"))
  #system("for i in *_R1_mm10.sam; do mv $i ${i%_R1_mm10.sam}.mm10.sam; done")
  #system("for i in *_R1_mm10.sam.log; do mv $i ${i%_R1_mm10.sam.log}.mm10.sam.log; done")
  system("mv *mm10.sam ../data/mapped_mm10")
  system("mv *sam.log ../data/stats_mapping")
}

print("bowtie2 mapping complete!")

