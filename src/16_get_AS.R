## This script is used to get AS, intronic.AS and intronic from ORF-removed transcripts
## Date: 18/02/2014
## Author: Zhipeng

library("Rsamtools")
library("GenomicRanges")

##########fun
rm(list = ls())
myfun = function(file){
      load("../DB/geneModel.RData")
      ###exonAS
      exonAS.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.AS.noBlastx.noORF.fa", sep = ""), format = "fasta", use.names = T)
      writeXStringSet(exonAS.fa, paste("../LincRNA/", file, ".assembled.exon.gene.exonAS.fa", sep = ""))
      load(paste("../Rdata/", file, ".assembled.annotation.grl.RData", sep = ""))
      exonAS.id = names(exonAS.fa)
      exonAS.grl = all.assembled.exon.gene.AS.grl[exonAS.id]
      exonAS.width = as.numeric(lapply(exonAS.grl, function(x) sum(width(x))))

      if(length(exonAS.fa) == 0){
      exonAS.neigh.df = NULL
      }else{
      exonAS.mRNA.gr = unlist(range(exonAS.grl))
      exonAS.mRNA.names = names(exonAS.grl)
      values(exonAS.mRNA.gr)$names = exonAS.mRNA.names 
      save(exonAS.grl, exonAS.mRNA.gr, file = paste("../Rdata/", file, ".exonAS.RData", sep = ""))

      exonAS.mRNA.ds.gr = exonAS.mRNA.gr
      strand(exonAS.mRNA.ds.gr) = rep("*", length(exonAS.mRNA.ds.gr))
      exonAS.ds.neigh = distanceToNearest(exonAS.mRNA.ds.gr, geneModel.mRNA.gr)
      exonAS.ds.neigh.nodup = exonAS.ds.neigh[!duplicated(queryHits(exonAS.ds.neigh))]
      gene.neigh.gr = geneModel.mRNA.gr[subjectHits(exonAS.ds.neigh.nodup)]
      gene.neigh.gr$distance = as.data.frame(exonAS.ds.neigh.nodup)$distance
      names(gene.neigh.gr) = names(exonAS.mRNA.ds.gr[queryHits(exonAS.ds.neigh.nodup)])
      gene.neigh.df = as.data.frame(gene.neigh.gr)
      exonAS.neigh.df = as.data.frame(exonAS.mRNA.gr)
      exonAS.neigh.df = data.frame(exonAS.neigh.df, gene.neigh.df[names(exonAS.mRNA.gr),])
      exonAS.neigh.df[,4] = exonAS.width
      }
      write.table(exonAS.neigh.df, file = paste("../LincRNA/", file, ".exonAS.neigh.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)

      ###intronAS
      intronAS.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.intronic.AS.noBlastx.noORF.fa", sep = ""), format = "fasta", use.names = T)
      writeXStringSet(intronAS.fa, paste("../LincRNA/", file, ".assembled.exon.gene.intronAS.fa", sep = ""))
      load(paste("../Rdata/", file, ".assembled.annotation.grl.RData", sep = ""))
      intronAS.id = names(intronAS.fa)
      intronAS.grl = all.assembled.exon.gene.intronic.AS.grl[intronAS.id]
      intronAS.width = as.numeric(lapply(intronAS.grl, function(x) sum(width(x))))

      if(length(intronAS.fa) == 0){
      intronAS.neigh.df = NULL
      }else{
      intronAS.mRNA.gr = unlist(range(intronAS.grl))
      intronAS.mRNA.names = names(intronAS.grl)
      values(intronAS.mRNA.gr)$names = intronAS.mRNA.names
      save(intronAS.grl, intronAS.mRNA.gr, file = paste("../Rdata/", file, ".intronAS.RData", sep = ""))

      intronAS.mRNA.ds.gr = intronAS.mRNA.gr
      strand(intronAS.mRNA.ds.gr) = rep("*", length(intronAS.mRNA.ds.gr))
      intronAS.ds.neigh = distanceToNearest(intronAS.mRNA.ds.gr, geneModel.mRNA.gr)
      intronAS.ds.neigh.nodup = intronAS.ds.neigh[!duplicated(queryHits(intronAS.ds.neigh))]
      gene.neigh.gr = geneModel.mRNA.gr[subjectHits(intronAS.ds.neigh.nodup)]
      gene.neigh.gr$distance = as.data.frame(intronAS.ds.neigh.nodup)$distance
      names(gene.neigh.gr) = names(intronAS.mRNA.ds.gr[queryHits(intronAS.ds.neigh.nodup)])
      gene.neigh.df = as.data.frame(gene.neigh.gr)
      intronAS.neigh.df = as.data.frame(intronAS.mRNA.gr)
      intronAS.neigh.df = data.frame(intronAS.neigh.df, gene.neigh.df[names(intronAS.mRNA.gr),])
      intronAS.neigh.df[,4] = intronAS.width
      }
      write.table(intronAS.neigh.df, file = paste("../LincRNA/", file, ".intronAS.neigh.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)

      ###intronic
      intronic.fa = readDNAStringSet(paste("../ORF/", file, ".assembled.exon.gene.intronic.noBlastx.noORF.fa", sep = ""), format = "fasta", use.names = T)
      writeXStringSet(intronic.fa, paste("../LincRNA/", file, ".assembled.exon.gene.intronic.fa", sep = ""))
      load(paste("../Rdata/", file, ".assembled.annotation.grl.RData", sep = ""))
      intronic.id = names(intronic.fa)
      intronic.grl = all.assembled.exon.gene.intronic.grl[intronic.id]
      intronic.width = as.numeric(lapply(intronic.grl, function(x) sum(width(x))))

      if(length(intronic.fa) == 0){
      intronic.neigh.df = NULL
      }else{
      intronic.mRNA.gr = unlist(range(intronic.grl))
      intronic.mRNA.names = names(intronic.grl)
      values(intronic.mRNA.gr)$names = intronic.mRNA.names
      save(intronic.grl, intronic.mRNA.gr, file = paste("../Rdata/", file, ".intronic.RData", sep = ""))

      intronic.mRNA.ds.gr = intronic.mRNA.gr
      strand(intronic.mRNA.ds.gr) = rep("*", length(intronic.mRNA.ds.gr))
      intronic.ds.neigh = distanceToNearest(intronic.mRNA.ds.gr, geneModel.mRNA.gr)
      intronic.ds.neigh.nodup = intronic.ds.neigh[!duplicated(queryHits(intronic.ds.neigh))]
      gene.neigh.gr = geneModel.mRNA.gr[subjectHits(intronic.ds.neigh.nodup)]
      gene.neigh.gr$distance = as.data.frame(intronic.ds.neigh.nodup)$distance
      names(gene.neigh.gr) = names(intronic.mRNA.ds.gr[queryHits(intronic.ds.neigh.nodup)])
      gene.neigh.df = as.data.frame(gene.neigh.gr)
      intronic.neigh.df = as.data.frame(intronic.mRNA.gr)
      intronic.neigh.df = data.frame(intronic.neigh.df, gene.neigh.df[names(intronic.mRNA.gr),])
      intronic.neigh.df[,4] = intronic.width
      }
      write.table(intronic.neigh.df, file = paste("../LincRNA/", file, ".intronic.neigh.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)

}

setwd("../../LincRNA")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
