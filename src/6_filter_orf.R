##This script is used to fiter transcripts with long ORFs
##Date: 18/02/2014
##Author: Zhipeng

library("Rsamtools")

#######function
rm(list = ls())

myfun = function(file){
      seqs.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.intergenic.noBlastx.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../ORF/", file, ".assembled.exon.gene.intergenic.noBlastx.withORF.id", sep = ""))$size == 0){
        orf = NULL
      } else {
        orf = read.delim(paste("../ORF/", file, ".assembled.exon.gene.intergenic.noBlastx.withORF.id", sep = ""), header = F)
      }
 
      seqs.noORF.fa = seqs.fa[!(names(seqs.fa) %in% unique(orf$V1))]
      writeXStringSet(seqs.noORF.fa, paste("../ORF/", file, ".assembled.exon.gene.intergenic.noBlastx.noORF.fa", sep = ""))
}

setwd("../../ORF")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
