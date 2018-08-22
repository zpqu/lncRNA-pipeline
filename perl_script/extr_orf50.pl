#! /usr/bin/perl
 
open(IN,$ARGV[0]);
open(IN1,$ARGV[1]);
open(OUT,">".$ARGV[2]);
open(OUT1,">".$ARGV[3]);
open(OUT2,">".$ARGV[4]);

$/=">";

while(<IN>){
    ~s/>//;
    next if ($_ eq "");
    @seq=split(/\n/,$_,2);
    $seq_head=$seq[0];
#    print $seq_head."\n";
    $seq[1]=~s/\n//g;
    $seq_length=length($seq[1]);

    if($seq_head=~/\|gb\|(\w+)/){
#    print $1."\n";
	$head=$1;
	$head_length{$head}=$seq_length+1;
	$full_id{$head}=$seq_head;
    }else{
	$seq_head=~/^([^\s]+)/;
	$seq_head_s=$1;
	$head_length{$seq_head_s}=$seq_length+1;
	$full_id{$seq_head_s}=$seq_head;
    }
#    print $head."\t".$head_length{$head}."\n";
}


while(<IN1>){
    ~s/>//;
    next if ($_ eq "");
    
    if($_=~/^([^\s]+)\_\d+\s\[(\d+)\s\-\s(\d+)/){
	$id=$1;
	$start=$2;
	$end=$3;
    }
#    print $id."\t".$start."\t".$end."\n";
    $query_length_abs=abs($start-$end);
    if ($query_length_abs >=300){
	print OUT ">".$_;
	print OUT1 $id."\t"."more than 100 AA"."\n";
	$hash_uni_id{$id}="Y";
    }
    if ($start < $end){
	$end_length=$head_length{$id}-$end;
	$query_length=$end - $start;
	print $_ if ($head_length{$id} == 0);
	$coverage=int($query_length*100/$head_length{$id});
	if  (($start < 4) || ($end_length < 3)){
	    print OUT ">".$_;
	    print OUT1 $id."\t".$start."\t".$end."\t".$head_length{$id}."\t".$coverage."\n";
	    $hash_uni_id{$id}="Y";
	}
    } else {
	$end_length=$head_length{$id}-$start;
	$query_length=$start - $end;
	$coverage=int($query_length*100/$head_length{$id});
	if (($end < 4 ) || ($end_length < 3)){
	    print OUT ">".$_;
	    print OUT1 $id."\t".$start."\t".$end."\t".$head_length{$id}."\t".$coverage."\n";
	    $hash_uni_id{$id}="Y";
	}
    }

}

foreach $orf_id (keys %hash_uni_id) {
    print OUT2 $full_id{$orf_id}."\n";
}


close IN;
close IN1;
close OUT;
close OUT1;
close OUT2;
