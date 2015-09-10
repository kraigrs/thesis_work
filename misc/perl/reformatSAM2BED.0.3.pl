#!/usr/bin/perl

########################################################################
# 
# 06/09/2010
#
# reformatSAM2BED.0.2.pl
#
# Purpose: makes a .bed file and lift the coordinates from one genome 
#          to another
# 
# Input: .sam alignment file from MOSAIK (or other alignment tool)
#
# Output: .bed file, and lifted .bed file
#
# Usage: perl reformatSAM2BED.0.2.pl <samfile> <lift[yes/no]> <chain>
#        <samfile> ==> .sam file
#        <lift>    ==> say whether or not to perform the liftover on the new .bed file     
#        <chain>   ==> liftover chain file, write NA if not applicable  
#
########################################################################

use strict;
use warnings;

main ();
sub main 
{
  my$samfile = $ARGV[0];
  my$lift = $ARGV[1];
  my$chain = $ARGV[2];
 
  my$sam; 
  my$read; 
  my$contig; 
  my$start; 
  my$end;
  my@fields;
  my@seq; 
  
  if($samfile =~ /(\S+)\.sam$/){$sam = $1;}

  open(IN,"$samfile") or die "Can't open $samfile for reading!";
  open(OUT,">$sam\.bed") or die "Error writing to $sam\.bed\n";

  HERE:while(<IN>)
  {
    chomp;
    if($_ =~ /^\@.+/){next HERE;}
    else
    {
      @fields = split(/\s+/,$_);  
      $read = $fields[0];
      $contig = $fields[2];
      $start = $fields[3] - 1; # mosaik may have a bug in outputting .bed files, in this case, must subtract 1
      @seq = split("",$fields[9]);
      $end = $fields[3] + scalar@seq; #fixed error, must subtract 1 ## may not actually have to subtract 1 ##
    
      # SAM files are always represented in the forward strand
      print OUT "$contig\t$start\t$end\t$read\n";
    } 
  }
  close IN;
  close OUT;

  if($lift eq "yes"){system "liftOver $sam\.bed $chain $sam\.bed.lifted $sam\.bed.lifted.un";}
}

