#!/usr/bin/perl

###################################################################################
#
# 06/09/2010
#
# convertChr2Gene.0.2.pl
#
# Purpose: query a gap file to see if mapped read falls in a gap and, if not, 
#          in which gene and exon it happens to fall in 
# 
# Input: .map file, list of gaps, list of constitutive exons
#
# Output: for every read, an indication of being in a gap, constitutive exon, or neither 
#
# Usage: perl convertChr2Gene.0.2.pl <bedfile> <gaps> <constExons> <output>  
#        <bedfile>    ==> the .bed file to be converted 
#        <gaps>       ==> a file containing the gap regions
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
  my%gaps;
  
  my$bedfile = $ARGV[0];  
  my$gaps_file = $ARGV[1];
  my$constEx_file = $ARGV[2]; 
  my$read_out = $ARGV[3];

  open (GAPS,"$gaps_file") or die "Can't open $gaps_file for reading!\n";
  while (<GAPS>) 
  {
    chomp;
    unless(/^track.+/)
    {
      @elements = split(/\s+/,$_); 	  
      $gaps{$elements[0]}{$elements[1]} = $elements[2];
      #gaps{chromosome}{start} = stop
    }
  }
  close GAPS;

  open (EXONS,"$constEx_file") or die "Can't open $constEx_file for reading!\n";
  while (<EXONS>) 
  {
    chomp;
    unless (/^track.+/)
    {
      @elements = split(/\s+/,$_); 
      $exons{$elements[0]}{$elements[1]} = [$elements[2], $elements[3]];
      #exons{chrm}{start} = [stop, gene]
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
    
    if ($line =~ /^track.+/){next LINE;} 
    else 
    {
      @elements = split(/\s+/,$line); 
      $chrom = $elements[0];
      $start = $elements[1];
      $stop = $elements[2];
      $read = $elements[3];
	  
      foreach $begin (keys %{$gaps{$chrom}}) 
      {
        $end = $gaps{$chrom}{$begin};
        #gaps{chromosome}{start} = stop
	    
        if ( ($stop >= $begin && $stop <= $end) || ($start >= $begin && $start <= $end) ) 
        {
          print OUT "gap_$chrom\t$begin\t$end\t$read\n";
          next LINE;
        }
      }
      foreach $begin (keys %{$exons{$chrom}}) 
      {
        #$exons{chrm}{start} = [stop, gene]
        #@elements = split(/,/,$coord);
        #$begin = $coord;
        $end = $exons{$chrom}{$begin}[0];
	    
        if ( ($stop >= $begin && $stop <= $end) || ($start >= $begin && $start <= $end ) ) 
        {
          print OUT "$exons{$chrom}{$begin}[1]\t$begin\t$end\t$read\n";
          next LINE;
        }
      }
      print OUT "noGap_noExon\t\*\t\*\t$read\n";
      ##This will only happen if neither the gap or exon requirements are met       
    } 
  }
  close BED;
  close OUT;
}
