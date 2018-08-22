## This script is used to get lincRNAs from ORF-removed transcripts
## Date: 18/02/2014
## Author: Zhipeng

library("Rsamtools")
library("GenomicRanges")

##########fun
rm(list = ls())
myfun = function(file){
      load("../DB/geneModel.RData")
      lincRNA.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.intergenic.noBlastx.noORF.fa", sep = ""), format = "fasta", use.names = T)
      writeXStringSet(lincRNA.fa, paste("../LincRNA/", file, ".assembled.exon.gene.lincRNA.fa", sep = ""))
      load(paste("../Rdata/", file, ".assembled.annotation.grl.RData", sep = ""))
      lincRNA.id = names(lincRNA.fa)
      lincRNA.grl = all.assembled.exon.gene.intergenic.grl[lincRNA.id]
      lincRNA.width = as.numeric(lapply(lincRNA.grl, function(x) sum(width(x))))
      
      if(length(lincRNA.fa) == 0){
      lincRNA.neigh.df = NULL
      }else{
      lincRNA.mRNA.names = names(lincRNA.grl)
      lincRNA.mRNA.gr = unlist(range(lincRNA.grl))
      values(lincRNA.mRNA.gr)$names = lincRNA.mRNA.names
      save(lincRNA.grl, lincRNA.mRNA.gr, file = paste("../Rdata/", file, ".lincRNA.RData", sep = ""))

      lincRNA.mRNA.ds.gr = lincRNA.mRNA.gr
      strand(lincRNA.mRNA.ds.gr) = rep("*", length(lincRNA.mRNA.ds.gr))
      lincRNA.ds.neigh = distanceToNearest(lincRNA.mRNA.ds.gr, geneModel.mRNA.gr)
      lincRNA.ds.neigh.nodup = lincRNA.ds.neigh[!duplicated(queryHits(lincRNA.ds.neigh))]
      gene.neigh.gr = geneModel.mRNA.gr[subjectHits(lincRNA.ds.neigh.nodup)]
      gene.neigh.gr$distance = as.data.frame(lincRNA.ds.neigh.nodup)$distance
      names(gene.neigh.gr) = names(lincRNA.mRNA.ds.gr[queryHits(lincRNA.ds.neigh.nodup)])
      gene.neigh.df = as.data.frame(gene.neigh.gr)
      lincRNA.neigh.df = as.data.frame(lincRNA.mRNA.gr)
      lincRNA.neigh.df = data.frame(lincRNA.neigh.df, gene.neigh.df[names(lincRNA.mRNA.gr),])
      lincRNA.neigh.df[,4] = lincRNA.width
      }
      write.table(lincRNA.neigh.df, file = paste("../LincRNA/", file, ".lincRNA.neigh.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
}

setwd("../../LincRNA")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
