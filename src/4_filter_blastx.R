## This script is used to fiter blastx results
## Date: 18/02/2014
## Author: Zhipeng

library("Rsamtools")

##########
rm(list = ls())

myfun = function(file){
      seqs.fa = readDNAStringSet(paste("../Output/", file, ".assembled.exon.gene.intergenic.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../blastx/", file, ".assembled.exon.gene.intergenic.blastx.1e-3.out", sep = ""))$size == 0){
        blastx = NULL
      } else {
        blastx = read.delim(paste("../blastx/", file, ".assembled.exon.gene.intergenic.blastx.1e-3.out", sep = ""), header = F)
      }
 
      seqs.noBlastx.fa = seqs.fa[!(names(seqs.fa) %in% unique(blastx$V1))]
      writeXStringSet(seqs.noBlastx.fa, paste("../ORF/", file, ".assembled.exon.gene.intergenic.noBlastx.fa", sep = ""))
}

setwd("../../ORF")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}


