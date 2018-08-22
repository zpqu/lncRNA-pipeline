##

cd ../../ORF

FILES=$(ls *noBlastx.fa)
for f in $FILES
do
    echo "Processing $f ..."
    filename=${f%.*}
    getorf -minsize 150 -sequence $f -outseq ${filename}.ORF50aa.fa ##add -reverse N if it's directional sequencing
    perl ../Script/perl_script/extr_orf50.pl $f ${filename}.ORF50aa.fa ${filename}.withORF.fa ${filename}.withORF.list ${filename}.withORF.id
done
