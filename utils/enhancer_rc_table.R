library("glue")
library("argparse")
set.cores = 20

# create parser object
parser = ArgumentParser()

parser$add_argument("-w", "--workdir", type = "character",
                    help = "name of working dir")
parser$add_argument("-b", "--bed", type = "character",
                    help = "bed file for subsetting")
parser$add_argument("-o", "--output", type = "character",
                    help = "output filename")
args = parser$parse_args()

# add working directory
workdir = args$workdir
bed = args$bed
output = args$output

setwd(paste0(workdir, "/mapped_mm10_20230510_EpiLC")) # path to alignment files

# index files
#system("for i in *_mm10.q10.fixmate.psort.bam; do samtools index $i; done")

#deeptools
print("multiBamSummary is in progress...")
system(paste0("multiBamSummary BED-file --BED ", bed, " --bamfiles *_mm10.q10.fixmate.psort.bam --outRawCounts ", output))