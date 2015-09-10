#! /usr/bin/perl

use strict;
use warnings;
#use 5.010;


our(@gene_info, @meta_info);
our($lastInfo1, $lastInfo2);

my $line;
my $lastMeta;
my $is_repeat = 0;
my $is_first = 1;

my $stackInfo1 = "No";
my $stackInfo2 = "No";

if(!open INPUT, "<$ARGV[0]"){
	die "Cannot open file $ARGV[0]";
}

open OUT, ">$ARGV[1]";

while(<INPUT>){
	chomp;
	$line = $_;
	@gene_info = split("\t",$line);
	@meta_info = split(/\@/,$gene_info[0]);
	if($is_first){
		$is_first = 0;
		$lastInfo1 = $gene_info[1];
		$lastInfo2 = $gene_info[2];
		$lastMeta = $meta_info[0];
		next;
	}
		
	if($meta_info[0] ne $lastMeta){
		print OUT $lastMeta."\t".$lastInfo1."\t".$lastInfo2."\n";
		$is_repeat = 0;
	}
	else{
		$is_repeat = 1;
		if($gene_info[1] ne "No" && $lastInfo1 eq "No"){
			$lastInfo1 = $gene_info[1];
		}
		if($gene_info[2] ne "No" && $lastInfo2 eq "No"){
			$lastInfo2 = $gene_info[2];
		}
		next;
	}
	$lastMeta = $meta_info[0];
	
	$lastInfo1 = $gene_info[1];
	$lastInfo2 = $gene_info[2];
}

print OUT $lastMeta."\t".$lastInfo1."\t".$lastInfo2."\n";

close INPUT;
close OUT;
