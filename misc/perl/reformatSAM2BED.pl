#!/usr/bin/perl

########################################################################
# 
# 04/01/2010
#
# reformatSAM2BED.pl
#
# Purpose: makes a .bed file and lift the coordinates from one genome 
#          to another
# 
# Input: a file from Bowtie (or other alignment tool)
#
# Output: .bed file
#
########################################################################

use strict;
use warnings;

my$mel1sam = $ARGV[0];
my$mel2sam = $ARGV[1];
my$sec1sam = $ARGV[2];
my$sec2sam = $ARGV[3];
my$melout = $ARGV[4];
my$secout = $ARGV[5];
my$secbed = $ARGV[6];
my$chain = $ARGV[7];

my$mel1line; my$mel2line; my$sec1line; my$sec2line;
my$mel1read; my$mel2read; my$sec1read; my$sec2read;
my$mel1contig; my$mel2contig; my$sec1contig; my$sec2contig;
my$mel1start; my$mel2start; my$sec1start; my$sec2start;
my$mel1end; my$mel2end; my$sec1end; my$sec2end;
my@mel1fields; my@mel2fields; my@sec1fields; my@sec2fields; 
my@mel1seq; my@mel2seq; my@sec1seq; my@sec2seq;
 
open(MEL1,"$mel1sam") or die "Can't open $mel1sam for reading!";
open(MEL2,"$mel2sam") or die "Can't open $mel2sam for reading!";

open(OUT,">$melout") or die "Error writing to $melout\n";

HERE:while(<MEL1>)
{
  $mel1line = $_;
  $mel2line = <MEL2>;
  chomp $mel1line;
  chomp $mel2line;

  if($mel1line =~ /^\@.+/){next HERE;}
  else
  {
    @mel1fields = split(/\s+/,$mel1line);  
    $mel1read = $mel1fields[0];
    $mel1contig = $mel1fields[2];
    $mel1start = $mel1fields[3];
    @mel1seq = split("",$mel1fields[9]);
    $mel1end = $mel1start + scalar@mel1seq - 1; #fixed error, must subtract 1

    @mel2fields = split(/\s+/,$mel2line);  
    $mel2read = $mel2fields[0];
    $mel2contig = $mel2fields[2];
    $mel2start = $mel2fields[3];
    @mel2seq = split("",$mel2fields[9]);
    $mel2end = $mel2start + scalar@mel2seq - 1; #fixed error, must subtract 1
    
    # SAM files are always represented in the forward strand
    print OUT "$mel1contig\t$mel1start\t$mel1end\t$mel1read\t\.\t\+\n"; 
    print OUT "$mel2contig\t$mel2start\t$mel2end\t$mel2read\t\.\t\+\n";
  }
}
close OUT;
close MEL1;
close MEL2;

open(SEC1,"$sec1sam") or die "Can't open $sec1sam for reading!";
open(SEC2,"$sec2sam") or die "Can't open $sec2sam for reading!";

open(OUT1,">$secout") or die "Error writing to $secout\n";
open(OUT2,">$secbed") or die "Error writing to $secbed\n";

HERE:while(<SEC1>)
{
  $sec1line = $_;
  $sec2line = <SEC2>;
  chomp $sec1line;
  chomp $sec2line;

  if($sec1line =~ /^\@.+/){next HERE;}
  else
  {
    @sec1fields = split(/\s+/,$sec1line);  
    $sec1read = $sec1fields[0];
    $sec1contig = $sec1fields[2];
    $sec1start = $sec1fields[3];
    @sec1seq = split("",$sec1fields[9]);
    $sec1end = $sec1start + scalar@sec1seq - 1;

    @sec2fields = split(/\s+/,$sec2line);  
    $sec2read = $sec2fields[0];
    $sec2contig = $sec2fields[2];
    $sec2start = $sec2fields[3];
    @sec2seq = split("",$sec2fields[9]);
    $sec2end = $sec2start + scalar@sec2seq - 1;
    
    # SAM files are always represented in the forward strand

    print OUT1 "$sec1contig\t$sec1start\t$sec1end\t$sec1read\t\.\t\+\n"; 
    print OUT1 "$sec2contig\t$sec2start\t$sec2end\t$sec2read\t\.\t\+\n";
    
    if($sec1contig ne "*")
    { 
      print OUT2 "$sec1contig\t$sec1start\t$sec1end\t$sec1read\t\.\t\+\n"; 
    }
    if($sec2contig ne "*")
    { 
      print OUT2 "$sec2contig\t$sec2start\t$sec2end\t$sec2read\t\.\t\+\n";
    } 
  }
}
close OUT1;
close OUT2;
close SEC1;
close SEC2;

system "./liftOver.MacOSX.ppc $secbed $chain $secbed\.lifted $secbed\.lifted.un";


