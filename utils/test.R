###:: mapping_mm10.sciTIP.R ::#####
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
parser$add_argument("-s", "--suffix", type = "character",
                    help = "add suffix string to the folders")
args = parser$parse_args()
# add folder containing FASTQ files
fastqs = args$input_dir
suffix = args$suffix

# bowtie2 index
btindex.path = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10_bowtie_index/mm10"
chromsize.path = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt"

system(paste0("mkdir -p ../data/mapped_mm10_", suffix))
system(paste0("mkdir -p ../data/stats_mapping_", suffix))