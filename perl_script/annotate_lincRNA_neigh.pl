#! /usr/bin/perl

open(IN,$ARGV[0]);
open(OUT,">".$ARGV[1]);

while(<IN>){
    chomp;
    next if ($_ eq "");
    next if ($_=~/^seqnames/);
    @tmp=split("\t",$_);
    $linc_strand=$tmp[4];
    $neigh_strand=$tmp[10];
    $linc_start=$tmp[2];
    $linc_end=$tmp[3];
    $neigh_start=$tmp[7];
    $neigh_end=$tmp[8];
    if($linc_strand eq "+" && $neigh_strand eq "+"){
	$type="Sense";
	if($linc_start < $neigh_start){
	    $end="5END";
	}else{
	    $end="3END";
	}
	print OUT $_."\t".$type."\t".$end."\n";
    }
    if($linc_strand eq "-" && $neigh_strand eq "-"){
	$type="Sense";
	if($linc_start < $neigh_start){
	    $end="3END";
	}else{
	    $end="5END";
	}
	print OUT $_."\t".$type."\t".$end."\n";
    }
    if($linc_strand eq "+" && $neigh_strand eq "-"){
	$type="AS";
	if($linc_start < $neigh_start){
	    $end="3END";
	}else{
	    $end="5END";
	}
	print OUT $_."\t".$type."\t".$end."\n";
    }
    if($linc_strand eq "-" && $neigh_strand eq "+"){
	$type="AS";
	if($linc_start < $neigh_start){
	    $end="5END";
	}else{
	    $end="3END";
	}
	print OUT $_."\t".$type."\t".$end."\n";
    }
    if($linc_strand eq "*"){
	$type="NA";
	if(($neigh_strand eq "+" && $linc_start < $neigh_start) || ($neigh_strand eq "-" && $linc_start > $neigh_start)){
	    $end="5END";
	}else{
	    $end="3END";
	}
	print OUT $_."\t".$type."\t".$end."\n";
    }
}

close IN;
close OUT;
