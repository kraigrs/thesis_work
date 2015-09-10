#!/usr/bin/perl

###################################################################################
#
# 04/14/2009
#
# cconvertChrom2Gene.pl
#
# Purpose: query a gap file to see if mapped read falls in a gap and, if not, 
#          in which gene and exon it happens to fall in 
# 
# Input: .map file, list of gaps, list of constitutive exons
#
# Output: for every read, an indication of being in a gap, constitutive exon, or neither 
#
# E.g. 
#
###################################################################################

use strict;
use warnings;

main ();

sub main 
{ 
  my @elements;
  my %exons;
  my %gaps;
  
  my $constEx_file = $ARGV[0];
  my $gaps_file = $ARGV[1];
  
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
  
  my $coord;
  my $line; 
  my $begin; 
  my $end;
  my $chrom; 
  my $start;
  my $stop;  
  my $read;
  
  my $bedfile = $ARGV[2];
  my $read_out = $ARGV[3]; 
  
  open(BED,"$bedfile") or die "Can't open $bedfile for reading!\n";
  open(OUT,">$read_out") or die "Error writing to $read_out!\n";
  
  LINE:while(<BED>)
  {
    chomp;
    $line = $_;
    
    if ($line =~ /^track.+/)
    {
      next LINE;
    } 
    else 
    {
      @elements = split(/\s+/,$line); 
      $chrom = $elements[0];
      $start = $elements[1];
      $stop = $elements[2];
      $read = $elements[3];
	  
      #check if read has a valid alignment
      if ($chrom eq "*" || $chrom eq "no_lift") 
      {
        print OUT "$chrom\t$start\t$stop\t$read\n";
        next LINE;
      } 
      else 
      {
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
  }
  close BED;
  close OUT;
}
