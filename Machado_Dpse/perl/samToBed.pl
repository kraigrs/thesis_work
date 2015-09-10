#!/usr/bin/perl

########################################################################
# 
# 02/21/2012
#
# samToBed.pl
#
# Purpose: makes a .bed file
# 
# Input: .sam alignment file 
#
# Output: .bed file
#
# Usage: perl samToBed.pl <samfile>
#        <samfile> ==> .sam file  
#
########################################################################

use strict;
use warnings;

main ();
sub main 
{
  my$samfile = $ARGV[0];

  my$sam; 
  my$read; 
  my$contig; 
  my$start; 
  my$end;
  my@fields;
  my@meta;
  my$seq; 
  
  if($samfile =~ /(\S+)\.sam$/){$sam = $1;}
  else{die "Error: please make sure you actually have a SAM file";}

  open(IN,"$samfile") or die "Can't open $samfile for reading!";
  open(OUT,">$sam\.bed") or die "Error writing to $sam\.bed\n";

  HERE:while(<IN>)
  {
    chomp;
    if($_ =~ /^\@.+/){next HERE;}
    else
    {
      @fields = split(/\t/,$_);
      unless($fields[2] eq "\*")
      {
        $read = $fields[0];
        $start = $fields[3] - 1; # SAM files are 1-based, so subtract 1 to get 0-based BED start
        $seq = $fields[9];
        $end = $start + length($seq);
        $contig = $fields[2];

        #@meta = split(/\|/,$fields[2]);
        #$contig = $meta[0];
    
        # SAM files are always represented in the forward strand
        print OUT "$contig\t$start\t$end\t$read\n";
        #print OUT "$contig\t$start\t$end\t$read\@$seq\@$start\@$end\n";
      }
    } 
  }
  close IN;
  close OUT;
}

