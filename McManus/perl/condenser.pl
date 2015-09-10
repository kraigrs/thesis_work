#!/usr/bin/perl

## first script condenses duplicates

$unique=0; 
while(<>) { 
	if (!($save{$_}++)) { 
		print $_; 
		$unique++;
		} 
		} 
		warn "\nChose $unique unique lines out of $. total lines.\n\n";

