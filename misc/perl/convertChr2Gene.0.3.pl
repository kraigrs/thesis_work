#!/usr/bin/perl

###################################################################################
#
# 06/09/2010
#
# convertChr2Gene.0.3.pl
#
# Purpose: query a list of constitutive exons to see where it falls 
# 
# Input: map file, list of constitutive exons
#
# Output: for every read an indication of being in a constitutive exon 
#
# Usage: perl convertChr2Gene.0.3.pl <bedfile> <constExons> <output>  
#        <bedfile>    ==> the .bed file to be converted 
#        <constExons> ==> a file containing a list of constitutive exons
#        <output>     ==> the output file to write converted reads to
#
###################################################################################

use strict;
use warnings;

main ();
sub main 
{ 
  my@elements;
  my%exons;
  
  my$bedfile = $ARGV[0];  
  my$constEx_file = $ARGV[1]; 
  my$read_out = $ARGV[2];

  open (EXONS,"$constEx_file") or die "Can't open $constEx_file for reading!\n";
  while (<EXONS>) 
  {
    chomp;
    unless (/^track.+/)
    {
      @elements = split(/\s+/,$_); 
      $exons{$elements[0]}{$elements[1]} = [$elements[2], $elements[3]];
      #exons{chrom}{start} = [stop, gene]
    }
  }
  close EXONS;

  my $coord;
  my $line; 
  my $begin; 
  my $end;
  my $chrom; 
  my $start;
  my $stop;  
  my $read;

  open(BED,"$bedfile") or die "Can't open $bedfile for reading!\n";
  open(OUT,">$read_out") or die "Error writing to $read_out!\n";
  
  LINE:while(<BED>)
  {
    chomp;
    $line = $_;
    
    if ($line =~ /^#track.+/){next LINE;} 
    else 
    {
      @elements = split(/\s+/,$line); 
      $chrom = $elements[0];
      $start = $elements[1];
      $stop = $elements[2];
      $read = $elements[3];
	  
      foreach $begin (keys %{$exons{$chrom}}) 
      {
        #$exons{chrom}{start} = [stop, gene]
        #@elements = split(/,/,$coord);
        #$begin = $coord;
        $end = $exons{$chrom}{$begin}[0];
	    
        if( ($stop >= $begin && $stop <= $end) || ($start >= $begin && $start <= $end) ) 
        {
          print OUT "$exons{$chrom}{$begin}[1]\t$begin\t$end\t$read\n";
          next LINE;
        }
      }       
    }
    print OUT "noExon\t\*\t\*\t$read\n"; 
  }
  close BED;
  close OUT;
}
