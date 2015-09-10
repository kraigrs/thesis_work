#!/usr/bin/perl

$col=1;
$limit=2;
$count=0;
while(<>) {
	s/\r?\n//;
	@F=split /\t/, $_;
	if ($F[$col] < $limit) {
		$count++;
		print "$_\n";
		}
	}
warn "\nChose $count lines out of $..\n\n";