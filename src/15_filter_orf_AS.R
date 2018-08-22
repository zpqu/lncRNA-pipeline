##This script is used to filter transcripts with long ORFs for AS, intronic.AS and intronic transcirpts
##Date: 12/05/2014
##Author: Zhipeng

library("Rsamtools")

#######function
rm(list = ls())

myfun = function(file){
      ###AS
      seqs.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.AS.noBlastx.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../ORF/", file, ".assembled.exon.gene.AS.noBlastx.withORF.id", sep = ""))$size == 0){
        orf = NULL
      } else {
        orf = read.delim(paste("../ORF/", file, ".assembled.exon.gene.AS.noBlastx.withORF.id", sep = ""), header = F)
      }

      seqs.noORF.fa = seqs.fa[!(names(seqs.fa) %in% unique(orf$V1))]
      writeXStringSet(seqs.noORF.fa, paste("../ORF/", file, ".assembled.exon.gene.AS.noBlastx.noORF.fa", sep = ""))

      ###AS_intronic
      seqs.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.intronic.AS.noBlastx.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../ORF/", file, ".assembled.exon.gene.intronic.AS.noBlastx.withORF.id", sep = ""))$size == 0){
        orf = NULL
      } else {
        orf = read.delim(paste("../ORF/", file, ".assembled.exon.gene.intronic.AS.noBlastx.withORF.id", sep = ""), header = F)
      }

      seqs.noORF.fa = seqs.fa[!(names(seqs.fa) %in% unique(orf$V1))]
      writeXStringSet(seqs.noORF.fa, paste("../ORF/", file, ".assembled.exon.gene.intronic.AS.noBlastx.noORF.fa", sep = ""))

      ###intronic
      seqs.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.intronic.noBlastx.fa", sep = ""), format = "fasta", use.names = T)
      if(file.info(paste("../ORF/", file, ".assembled.exon.gene.intronic.noBlastx.withORF.id", sep = ""))$size == 0){
        orf = NULL
      } else {
        orf = read.delim(paste("../ORF/", file, ".assembled.exon.gene.intronic.noBlastx.withORF.id", sep = ""), header = F)
      }

      seqs.noORF.fa = seqs.fa[!(names(seqs.fa) %in% unique(orf$V1))]
      writeXStringSet(seqs.noORF.fa, paste("../ORF/", file, ".assembled.exon.gene.intronic.noBlastx.noORF.fa", sep = ""))
}

setwd("../../ORF")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
