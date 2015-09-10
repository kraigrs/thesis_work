#!/usr/bin/perl
use strict;
use FileHandle;

my $line;
my $seq;
my $qual;
my $trim;
my $length;

open(READS, "<$ARGV[0]") or die("can't open reads file");  ## This is the varFiltered snp file (final snp file from samtools filtering)
while ($line =<READS>) {
    printf("$line");
    $line = <READS>;
    chomp($line);
    $seq = "$line";
    $trim = substr($seq, 0, $ARGV[1]);
	$length = length($trim);
    printf("$trim\n");
    $line = <READS>;
    chomp($line);
    printf("$line\n");
	<READS>;
    $qual = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    $trim = substr($qual, 0, $length);
    printf("$trim\n");
}
