#!/usr/bin/perl

$col=3;  
while (<>) { 
	s/\r?\n//; 
	@F = split /\t/, $_;
	 $val = $F[$col]; 
	 if (! exists $count{$val}) {
	 	 push @order, $val;
	 } 
	 	 $count{$val}++;
} 
foreach $val (@order) {
	print "$val\t$count{$val}\n";
	} 
warn "\nPrinted number of occurrences for ", scalar(@order), " values in $. lines.\n\n";