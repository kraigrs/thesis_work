#!/usr/bin/perl

use strict;
use warnings;

main ();
sub main 
{ 
  my@elements;
   
  my$constExons = $ARGV[0];
  my$outfile = $ARGV[1]; 

  open(OUT,">$outfile") or die "Error writing to $outfile!\n";

  open (EXONS,"$constExons") or die "Can't open $constExons for reading!\n";
  while (<EXONS>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      print OUT "$elements[0]\t$elements[1]\t$elements[2]\t$elements[3]\_$elements[1]\,$elements[2]\t$elements[4]\t$elements[5]\n";
    }
  }
  close EXONS;
  close OUT;
}
