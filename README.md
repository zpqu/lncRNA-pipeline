lncRNA identification pipeline

22/08/2018
Zhipeng

This pipeline is used to detect lncRNAs/lincRNAs from NGS data (RNA-seq) using bash, R and perl scripts.

Input files
Genomic coordinates (gtf) of transcripts (e.g. de novo assembled transcripts from cufflinks)

Output files
1) Identified genomic coordiante file (gtf) of lncRNAs, including lincRNAs (long intergenic non-coding RNAs), exonAS (exonic antisense transcripts), intronAS (intronic antisense transcripts), intronic (intronic non-coding RNAs).
2) Sequences of lncRNAs in fasta format.
3) Neighbour gene information for lincRNAs.
4) Summary table. 
