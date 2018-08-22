##This script is for retrieve lincRNAs from filtered transcripts
##Date:18/02/2014
##Author: Zhipeng 

cd ../../LincRNA
FILES=$(ls *.neigh.txt)
for f in $FILES
do
    filename=${f%.*}
    perl ../Script/perl_script/annotate_lincRNA_neigh.pl $f ${filename}.annotated.txt 
done
