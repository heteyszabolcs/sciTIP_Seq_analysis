# sciTIP-Seq analysis
Scripts for sciTIP-Seq data processing

Data processing scripts for H3.3 mESC sciTIP-Seq experiment.
The majority of the code is adapted from the original [Bartlett et al. paper](https://rupress.org/jcb/article/220/12/e202103078/212821/High-throughput-single-cell-epigenomic-profiling)

Additional codes: [Bartlett et al. github repo](https://github.com/dbart1807/TIP-seq)

### Steps: 
1. relocalize FASTQ files of each 96-well plate row. 
2. bowtie2 alignment to mm10 genome
3. BAM convertion and deduplications (samtools)
4. T7 deduplications
5. computing read count matrix (bedtools) 

Important folders at Uppmax:

  * Linxuan folder:
  _/proj/snic2020-6-3/LINXUAN/sciTIP_sc_analysis/_

  * Working directory:
  _/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils_
