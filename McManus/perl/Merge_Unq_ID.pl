#!/usr/bin/perl

$col1=0;
$col2=3;
($f1,$f2)=@ARGV;
open(F2,$f2);
while (<F2>) {
	 s/\r?\n//; 
	 @F=split /\t/, $_; 
	 $line2{$F[$col2]} .= "$_\n";
	 }
$count2 = $.;
open(F1,$f1);
while (<F1>) {
	s/\r?\n//;
	@F=split /\t/, $_;
	$x = $line2{$F[$col1]};
	if ($x) {
	  $num_changes = ($x =~ s/^/$_\t/gm);
	 	print $x; $merged += $num_changes;
	 	}
	 }
	 
warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n";