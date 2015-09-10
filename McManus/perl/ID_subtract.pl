#!/usr/bin/perl

($file1, $file2) = @ARGV; 
$printed = 0; 
open F2, $file2; 
while (<F2>) {
	 $h2{$_}++ 
}; 
$count2 = $.; 
open F1, $file1; 
while (<F1>) {
	 if (! $h2{$_}) {
	 print $_;
	 $printed++;
	 }
}
$count1 = $.; 
warn "\nRead $count1 lines from $file1 and $count2 lines from $file2.\nPrinted $printed lines found in $file1 but not in $file2\n\n";