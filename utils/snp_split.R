####################################
### Allele specific alignment of reads
### run SNPsplit on algined files
####################################

library("glue")
library("argparse")

# create parser object
parser = ArgumentParser()

parser$add_argument("-i", "--input_dir", type = "character",
                    help = "input directory with the processed BAM files")
parser$add_argument("-w", "--workdir", type = "character",
                    help = "name of working dir")
parser$add_argument("-o", "--output_dir", type = "character",
                    help = "output directory for SNPsplit outputs")
args = parser$parse_args()

# add folder containing deindexed FASTQ files
bams = args$input_dir
# add working directory
workdir = args$workdir
# add output directory
output = args$output_dir

# run SNPsplit
# make output directory
system(paste0("mkdir ", output))
setwd(workdir)

snpsplit = function(processed_bam, vcf_file = "CAST_EiJ_chr.vcf") {
  snpsplit = glue("SNPsplit --snp_file {vcf_file} {processed_bam}")
  system(snpsplit)
  system(paste0("mv *allele_flagged.bam ", output))
  system(paste0("mv *genome1.bam ", output))
  system(paste0("mv *genome2.bam ", output))
  system(paste0("mv *unassigned.bam ", output))
  system(paste0("mv *.SNPsplit_* ", output))
}

## loop over the finalized BAM files
# finding unprocessed BAMs
id = unlist(unname(sapply(list.files(bams, full.names = FALSE, pattern = "*.sam.final.bam"),
  strsplit, split = ".sam.final.bam")))
ready = unlist(unname(sapply(list.files(output, full.names = FALSE, pattern = "*.sam.final.genome1.bam"),
                             strsplit, split = "*.sam.final.genome1.bam")))
not_ready = setdiff(id, ready)

processed_bams = list.files(bams, full.names = TRUE, pattern = "*.sam.final.bam")
not_ready = sapply(not_ready, function(x) { grep(x, processed_bams) } )

# keeping BAMs with higher read counts
samtools_count = function(bam) {
  read_count = glue("samtools view -c {bam}")
  read_count = as.numeric(system(read_count, intern = TRUE))
  if (read_count > 1000) {
    return(bam)
  }
}
abundant_processed_bams = sapply(processed_bams[not_ready], samtools_count)
print("Kept only bams with above 1000 reads")
abundant_processed_bams[sapply(abundant_processed_bams, is.null)] = NULL

print(tail(abundant_processed_bams))
print(length(abundant_processed_bams))
sapply(abundant_processed_bams, snpsplit)



