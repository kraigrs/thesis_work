#!/usr/bin/perl

@cols=(0);  
while(<>) { 
	s/\r?\n//; 
	@F=split /\t/, $_; 
	print join("\t", @F[@cols]), "\n" 
} 

warn "\nChose columns ", join(", ", @cols), " for $. lines\n\n";