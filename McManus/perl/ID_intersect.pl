#!/usr/bin/perl

($file1, $file2) = @ARGV;
open F2, $file2 or die $!;
while (<F2>) {
	$h2{$_}++ 
}; 
open F1, $file1 or die;
$total=$.;
$printed=0;
while (<F1>) {
	$total++;
	if ($h2{$_}) {
		print $_;
		$h2{$_} = "";
		$printed++;
		}
	}
warn "\n\nRead $total lines.\nTook intersection and then removed duplicates, yielding $printed lines.\n";