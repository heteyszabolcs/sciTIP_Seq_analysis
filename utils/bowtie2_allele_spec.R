####################################
### Allele specific alignment of reads
### mapping to N-masked reference genome (coming from SNPsplit package)
####################################

library("travis") # coming from https://github.com/dvera/travis - see dependencies!
library("glue")
library("argparse")
set.cores = 20

# create parser object
parser = ArgumentParser()

parser$add_argument("-i", "--input_dir", type = "character",
                    help = "input directory with the fastq files")
parser$add_argument("-w", "--workdir", type = "character",
                    help = "name of working dir")
parser$add_argument("-o", "--output_dir", type = "character",
                    help = "output directory for the alignment files")
parser$add_argument("-b", "--bowtie_index", type = "character",
                    help = "path to folder containing bowtie index")
args = parser$parse_args()

# add folder containing deindexed FASTQ files
fastqs = args$input_dir
# add working directory
workdir = args$workdir
# add output directory
output = args$output_dir
# add bowtie index path
index = args$bowtie_index

# bowtie2 index
btindex.path = index

# # make output directory
system(paste0("mkdir ", output))

# # list fastqs
fq = unique(strsplit(
  list.files(
    path = fastqs,
    pattern = "R1_001.fastq.gz",
    full.names = TRUE
  ),
  '_R1_001.fastq.gz'
))

# chunk id
fq_id = unlist(unique(strsplit(
  list.files(
    path = fastqs,
    pattern = "R1_001.fastq.gz",
    full.names = FALSE
  ),
  '_R1_001.fastq.gz'
)))

cell_number = length(fq)

for (i in 1:length(fq)) {

  sams = unlist(strsplit(
    list.files(
      output,
      full.names = FALSE
    ),
    ".sam"
  ))

  # if fastq has been aligned already then skip
  print(glue("cell #{as.character(i)}"))
  if (glue("{fq_id[i]}") %in% sams) {
    next
  }

  print(paste0(as.character(i), " / ", as.character(cell_number)))

  bowtie2 = glue(
    "bowtie2 -x {btindex.path} -1 {fq[i]}_R1_001.fastq.gz -2 {fq[i]}_R2_001.fastq.gz --output {fq[i]}.sam --very-sensitive --end-to-end -N1 -q"
  )
  system(bowtie2)
  system(paste0("mv ",
                fastqs, "/", "*.sam ",
                output))
  system(paste0("mv ", fastqs, "/", "*sam.log ", output))
}

print("allele spec bowtie2 mapping complete!")

### Run after mapping script

####################################################################
######:: Remove duplicates -- samtools fixmate and markdup ::#######
####################################################################

setwd(output)
system(paste0("find . -type f -size 0 -delete")) # remove files with 0 byte
set.cores = 20

sams = unlist(strsplit(list.files(getwd(), 
                                  full.names = FALSE,
                                  pattern = ".sam"), ".sam"))

bams = unlist(strsplit(list.files(getwd(), 
                                  full.names = FALSE,
                                  pattern = ".final.bam"), ".sam.final.bam"))

sams = unname(sapply(setdiff(sams, bams), function(x) paste0(x, ".sam")))
print(glue("{length(sams)} are not converted to bam"))

##:: filter out <q10 reads and convert to name sorted BAM
samtools = function(sam) {
  
  print(sam)
  
  nsortBam = paste0("samtools view -q 10 -@ ", set.cores, " -u ", sam, " | samtools sort -n -o ", sam, ".q10.nsort.bam"
  )
  print("Working on bam converting and name sorting")
  system(nsortBam)
  
  ##:: run fixmate
  fixBam = paste0("samtools fixmate -@ ", set.cores, " -m ", sam, ".q10.nsort.bam ", sam, ".q10.nsort.fixmate.bam")
  print("Working on fixmate")
  system(fixBam)
  
  ##:: sort by coordinate
  psortBam = paste0("samtools sort -@ ", set.cores, " -o ", sam, ".q10.nsort.fixmate.psort.bam ", sam, ".q10.nsort.fixmate.bam")
  print("Working on coordinate sorting")
  system(psortBam)
  
  ##:: run markdup
  mkdupBams = paste0("samtools markdup -r -@ ", set.cores, " ", sam, ".q10.nsort.fixmate.psort.bam ", sam, ".q10.nsort.fixmate.psort.markdup.bam"
  )
  print("Working on duplicate marking")
  system(mkdupBams)
  
  # name sort (for bamToBed conversion)
  nsort_mkdupBams = paste0("samtools sort -n -@ ", set.cores," -o ", sam, ".final.bam " , sam, ".q10.nsort.fixmate.psort.markdup.bam"
  )
  print("Working on name sorting")
  system(nsort_mkdupBams)
}

# prepare bam files
sapply(sams, samtools)

print("samtools fixmate & markdup complete!")
