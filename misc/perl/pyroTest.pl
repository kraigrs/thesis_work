#!/usr/bin/perl

use strict;
use warnings;

main ();
sub main 
{
  my$ref1 = $ARGV[0];

  my$chr; my$segment; my$seqFound = 0;
  my%ref1; 

  open(REF1,"$ref1") or die "Can't open $ref1 for reading!\n";
  while(<REF1>) 
  {
    chomp;
    if(/^>(\S+)/)
    {
      if($seqFound == 1)
      {
        $ref1{$chr} = $segment;
        $seqFound = 0;
        undef $segment;
      }
      $chr = $1;
    }
    else
    {
      $seqFound = 1;
      $segment = $segment . $_;
      $ref1{$chr} = $segment;
    }
  }
  close REF1;
  undef $segment; $seqFound = 0;

  foreach(keys %ref1){print "$_\t$ref1{$_}\n";}
}
