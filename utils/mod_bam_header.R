library("GenomicAlignments")
library("rtracklayer")

param <- ScanBamParam(what=c("qname","flag","rname", "strand", "pos", "qwidth", "mapq", "cigar","mrnm","mpos","isize","seq","qual"),tag=c("NM","MD","MC","AS","XS"))

gr.aln <- readGAlignments("../data/cluster0.sorted.bam", 
                          param=param, use.names=TRUE)



add_id = function(x, id = "cluster0") {
  mod = paste0(id, ":", x)
  return(mod)
}

new_qname = sapply(mcols(gr.aln)$qname, add_id)
new_qname = unname(new_qname)
head(new_qname)

mcols(gr.aln)$qname = new_qname
export(gr.aln, "../data/cluster0_mod.sorted.bam")
