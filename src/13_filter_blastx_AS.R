## This script is used to filter blastx results for AS, intronic.AS and intronic transcripts
## Date: 12/05/2014
## Author: Zhipeng

library("Rsamtools")

##########
rm(list = ls())

myfun = function(file){
      seqs.fa = readDNAStringSet(paste("../Output/", file, ".assembled.exon.gene.AS.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../blastx/", file, ".assembled.exon.gene.AS.blastx.1e-3.out", sep = ""))$size == 0){
        blastx = NULL
      } else {
	blastx = read.delim(paste("../blastx/", file, ".assembled.exon.gene.AS.blastx.1e-3.out", sep = ""), header = F)
      }
      seqs.noBlastx.fa = seqs.fa[!(names(seqs.fa) %in% unique(blastx$V1))]
      writeXStringSet(seqs.noBlastx.fa, paste("../ORF/", file, ".assembled.exon.gene.AS.noBlastx.fa", sep = ""))

      seqs.fa = readDNAStringSet(paste("../Output/", file, ".assembled.exon.gene.intronic.AS.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../blastx/", file, ".assembled.exon.gene.intronic.AS.blastx.1e-3.out", sep = ""))$size == 0){
        blastx = NULL
      } else {
	blastx = read.delim(paste("../blastx/", file, ".assembled.exon.gene.intronic.AS.blastx.1e-3.out", sep = ""), header = F)
      }
      seqs.noBlastx.fa = seqs.fa[!(names(seqs.fa) %in% unique(blastx$V1))]
      writeXStringSet(seqs.noBlastx.fa, paste("../ORF/", file, ".assembled.exon.gene.intronic.AS.noBlastx.fa", sep = ""))

      seqs.fa = readDNAStringSet(paste("../Output/", file, ".assembled.exon.gene.intronic.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../blastx/", file, ".assembled.exon.gene.intronic.blastx.1e-3.out", sep = ""))$size == 0){
        blastx = NULL
      } else {
	blastx = read.delim(paste("../blastx/", file, ".assembled.exon.gene.intronic.blastx.1e-3.out", sep = ""), header = F)
      }
      seqs.noBlastx.fa = seqs.fa[!(names(seqs.fa) %in% unique(blastx$V1))]
      writeXStringSet(seqs.noBlastx.fa, paste("../ORF/", file, ".assembled.exon.gene.intronic.noBlastx.fa", sep = ""))
}

setwd("../../ORF")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
