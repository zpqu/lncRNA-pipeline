## This script is used to run blastx for intergenic transcripts against Swiss-port
## Date: 17/02/2014
## Author: Zhipeng

cd ../../blastx

echo -ne "The process of blastx against Uniprot_Sprot starts at: "
date

FILES=$(ls ../GTF/)
for f in $FILES
do
    echo -ne "Processing $f ..."
    date
    filename=${f%.*}
    blastx -outfmt 6 -db ../DB/uniprot_sprot.fasta -query ../Output/${filename}.assembled.exon.gene.intergenic.fa -out ${filename}.assembled.exon.gene.intergenic.blastx.1e-3.out -evalue 1e-3 -num_threads 8 -num_alignments 1
done

echo -ne "The process of blastx against Uniprot_Sprot ends at: "
date