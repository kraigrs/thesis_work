#!/usr/bin/perl

########################################################################
# 
# 03/27/2013
#
# generate_reads.pl
#
# Purpose: generate tiled sequence reads directly from FASTA file
#
########################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$length = $ARGV[1];
my$out = $ARGV[2];
 
my$first = 1; my$seq = "";
my$chr; my$pos; my$code; my$begin; my$end; my$temp; my$i; my$size;
my@elements; my@parts;
my%genome; my%rev;

print "Inputting $ref\n";

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #if(/^\>(\S.chr\S.)$/)
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome{$chr} = $seq;} # the first time, don't save sequence in hash
    $chr = $1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $chr\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

foreach $chr (keys %genome)
{
  $temp = &revComp($genome{$chr});
  $rev{$chr} = $temp;
}

print "Generating reads\n\n";

open(OUT," > $out") or die "\nError writing to $out\n";

foreach $chr (keys %genome)
{ 
  # positive orientation

  $size = length($genome{$chr});
  for($i=0;$i<=$size-$length;$i++)
  {
    $temp = substr $genome{$chr}, $i, $length;
    $begin = $i; $end = $i+$length; #BED format

    print OUT "\>$chr\|$begin\_$end\_pos\n";
    print OUT "$temp\n";
  }

  # negative orientation

  $size = length($rev{$chr});
  for($i=0;$i<=$size-$length;$i++)
  {
    $temp = substr $rev{$chr}, $i, $length;
    $begin = $size-$i-$length; $end = $size-$i;

    print OUT "\>$chr\|$begin\_$end\_neg\n";
    print OUT "$temp\n";
  }
}
close OUT;


sub revComp
{
  my$val = shift;
  my$revComp = reverse $val;
  $revComp =~ tr/ACGT/TGCA/;
  return $revComp;
}
