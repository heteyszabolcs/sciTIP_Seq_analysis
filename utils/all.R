####################################
### run after deindexing script
####################################

library("travis")
library("glue")
library("argparse")
set.cores = 20

# create parser object
parser = ArgumentParser()

parser$add_argument("-i", "--input_dir", type = "character",
                    help = "input directory with the fastq files")
parser$add_argument("-w", "--workdir", type = "character",
                    help = "name of working dir")
args = parser$parse_args()

# add folder containing deindexed FASTQ files
fastqs = args$input_dir
# add working directory
workdir = args$workdir

# bowtie2 index
btindex.path = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10_bowtie_index/mm10"
chromsize.path = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt"

system(paste0("mkdir -p ", workdir))
system(paste0("mkdir -p ", workdir, "/mapped_mm10"))
system(paste0("mkdir -p ", workdir, "/stats_mapping"))

fq = unique(strsplit(
  list.files(
    path = fastqs,
    pattern = "R1_001.fastq.gz",
    full.names = TRUE
  ),
  '_R1_001.fastq.gz'
))
print(length(fq)) # number of cells

for (i in 1:length(fq)) {
  Sam = bowtie2(
    read1files = paste0(fq[i], "_R1_001.fastq.gz"),
    read2files = paste0(fq[i], "_R2_001.fastq.gz"),
    indexfile = btindex.path,
    alignMode = "very-sensitive-local",
    appendIndexToName = TRUE,
    unaligned = FALSE,
    reorder = TRUE,
    threads = set.cores,
    minInsertSize = 10,
    maxInsertSize = 700
  )
  #system(glue("bowtie2 --threads 8 --very-sensitive-local --minins 10 --maxins 700 -x {btindex.path} -1 {fq[i]}_R1_001.fastq.gz -2 {fq[i]}_R2_001.fastq.gz -S {fq[i]}_mm10.sam"))
  #system("for i in *_R1_mm10.sam; do mv $i ${i%_R1_mm10.sam}.mm10.sam; done")
  #system("for i in *_R1_mm10.sam.log; do mv $i ${i%_R1_mm10.sam.log}.mm10.sam.log; done")
  system(paste0("mv *mm10.sam ", workdir, "/mapped_mm10"))
  system(paste0("mv *sam.log ", workdir, "/stats_mapping"))
}

print("bowtie2 mapping complete!")

######:: removeDups.R ::##

### Run after mapping script

####################################################################
######:: Remove duplicates -- samtools fixmate and markdup ::#######
####################################################################

library("glue")
library("travis")

setwd(paste0(workdir, "/mapped_mm10"))
set.cores = 20

##:: filter out <q10 reads and convert to name sorted BAM
nsortBam = paste0(
  "for i in *mm10.sam; do samtools view -q 10 -@ ",
  set.cores,
  " -u $i | samtools sort -n -o ${i%sam}q10.nsort.bam; done"
)
head(nsortBam)  ## delete
system(nsortBam)

### enable only if bam conversion successful!!
#system("rm *.sam")

##:: run fixmate
fixBam = paste0(
  "for i in *q10.nsort.bam; do samtools fixmate -@ ",
  set.cores,
  " -m $i ${i%nsort.bam}fixmate.bam; done"
)
system(fixBam)

##:: sort by coordinate
psortBam = paste0("for i in *fixmate.bam; do samtools sort -@ ",
                  set.cores,
                  " -o ${i%bam}psort.bam $i; done")
system(psortBam)

##:: run markdup
mkdupBams = paste0(
  "for i in *fixmate.psort.bam; do samtools markdup -r -@ ",
  set.cores,
  " $i ${i%mate.psort.bam}mkdup.bam; done"
)
system(mkdupBams)

# name sort (for bamToBed conversion)
nsort_mkdupBams = paste0(
  "for i in *.q10.fixmkdup.bam; do samtools sort -n -@ ",
  set.cores,
  " -o ${i%bam}nsort.bam $i; done"
)
system(nsort_mkdupBams)

print("samtools fixmate & markdup complete!")

##:: bamToBed conversion
nsort_mkdupBams = list.files(pattern = "*mkdup.nsort.bam")
Beds = bamToBed(
  nsort_mkdupBams,
  paired = TRUE,
  threads = set.cores,
  sortThreads = set.cores
)  # requires nsort bam

print("BAM to BED conversion complete!")

system(paste0("mkdir -p ", workdir, "/beds_rmdup"))
system(paste0("for i in *nsort.bed; do mv $i ../beds_rmdup", "/${i%nsort.bed}bed; done"))

print("samtools deduplication complete!")

### CLEANUP: enable only if fixmate.mkdup conversion successful!! to clear disk space if necessary
# system("rm *.q10.fixmate*")
# system("rm *.q10.fixmkdup.bam")

####:: removeT7dups.R ::####
#############################

####################################################################
######:: Remove T7 duplicates -- based on read1 start position ::#######
####################################################################

library("glue")
library("travis")

setwd(paste0(workdir, "/beds_rmdup"))

system("find . -type f -size 0 -delete") ## remove cells with 0 bytes to avoid error below
allBEDs = list.files(pattern = ".bed")
for (f in allBEDs) {
  tip = read.delim(f, header = F)
  tip2 = tip[, c(1, 2, 3)] ### will need to change this is using bedpe input --> tip2=tip[,c(1,2,6)]
  names(tip2) = c("CHR", "StartR1", "EndR2")
  tip2$ID1 = paste(tip2$CHR, tip2$StartR1, sep = "_")
  tip2$ID2 = paste(tip2$V1, tip2$EndR2, sep = "_")
  tip2$index = 1:nrow(tip2)
  undups = tip2[!duplicated(tip2$ID1), c("index")]
  tip3 = tip[undups, ]
  write.table(
    tip3,
    file = paste0(f, "_rmT7dup.bed"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    col.names = FALSE
  )
  cat(paste("Iteration", f, "was finished.\n"))
}

system("for i in *bed_rmT7dup.bed; do mv $i ${i%bed_rmT7dup.bed}rmT7dup.bed; done")
system(paste0("mkdir -p ", workdir, "/beds_rmT7dup"))
system(paste0("mv *.rmT7dup.bed ../beds_rmT7dup/"))

print("remove T7 duplications complete!")

####:: matrix.coverage.R ::####

##################################################
####:: calculate coverage and build matrix ::#####
##################################################

library(travis)

set.cores = 20
build = "mm10"
bin_size = 5000
step_size = 5000
scalar = 1

setwd(paste0(workdir, "/beds_rmdup")) # path to beds
allBEDs <- list.files(pattern = "fixmkdup.bed")

windows <-
  bedtoolsMakeWindows(
    bedfiles = chromsize.path,
    windowsize = bin_size,
    stepsize = step_size,
    threads = set.cores,
    genome = TRUE,
    outnames = paste0("~/bin/", build, "_win", bin_size, ".step", step_size, ".bed")
  )

system("for i in *bed; do bedSort $i $i; done")
cov = paste0("for i in *.bed; do bedtools coverage -counts -a ",
             windows,
             " -b $i > ${i%bed}win5k.bg; done")
system(cov)

###:: bigwig conversion (optional, for IGV)  # may need to run in command line...
bg2bw = paste0("for j in *.bg; do bedGraphToBigWig $j ",
               chromsize.path,
               " ${j%bg}bw; done")
system(bg2bw)
system(paste0("mkdir -p ", workdir, "/coverage_win5k"))
system(paste0("mv *.win5k.bg ", workdir, "/coverage_win5k"))
# system("mkdir -p ../coverage_win5k/win5k_bigwigs")
# system("mv *.win5k.bw ../coverage_win5k/win5k_bigwigs")

###########################################################
########--- build single-cell coverage matrix --- #########
###########################################################
set.cores = 20
setwd(paste0(workdir, "/coverage_win5k"))
allBgs <- list.files(pattern = "win5k.bg")
coords = read.tsv(allBgs[1])

chunk_length = 100
chunks = split(allBgs, ceiling(seq_along(allBgs) / chunk_length))
for (i in 1:length(chunks)) {
  print(paste0("chunk ", as.character(i), "/", as.character(chunk_length)))
  a = bgRead(unname(unlist(chunks[i])), threads = set.cores)
  fdata = cbind(coords[, 1:3], a)
  nam = names(fdata)
  nam = nam[4:length(nam)]
  nam = unlist(strsplit(nam, ".mm10.q10"))[c(TRUE, FALSE)]
  nam = c("CHR", "START", "END", nam)
  names(fdata) = nam
  write.table(
    fdata,
    file = gzfile(
      paste0("sciTIP_counts_matrix_chunk", as.character(i), ".txt.gz")
    ),
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
  
}

##:: cleanup ::##
system("rm *.win5k.bg") # optional: delete bedgraphs, do not need the after creating matrix and making bigwigs.
