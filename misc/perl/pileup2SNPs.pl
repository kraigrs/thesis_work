#!/usr/bin/perl

use strict;
use warnings;

my$pileup = $ARGV[0];

my$chr; my$pos; my$covseqs; my$call; my$match; my$As; my$Cs; my$Gs; my$Ts; my$mismatch; my$ref;
my@elements; my@bases; my@quals;

print "chrom\tpos\tref\tcov\tmatch\tmismatch\tA\tC\tG\tT\n";

open(PILE,"$pileup") or die "\nError opening $pileup\n";
while(<PILE>)
{
  chomp;
  unless(/^#/)
  {
    @elements = split(/\s+/,$_);
    $chr = $elements[0];
    $pos = $elements[1];
    $ref = $elements[2];
    $covseqs = $elements[3];
    $call = $elements[4];
  
    $match = 0;
    $As = 0;
    $Cs = 0;
    $Gs = 0;
    $Ts = 0; 
    $mismatch = 0;

    @bases = split("",$call);
    foreach(@bases)
    {
      if(/[\.,]/){$match+=1;}
      elsif(/[Aa]/){$As+=1;}
      elsif(/[Cc]/){$Cs+=1;}
      elsif(/[Gg]/){$Gs+=1;}
      elsif(/[Tt]/){$Ts+=1;}
    }
  
    $mismatch = $As + $Cs + $Gs + $Ts;

    if($ref eq "A"){$As += $match;}
    elsif($ref eq "C"){$Cs += $match;}
    elsif($ref eq "G"){$Gs += $match;}
    elsif($ref eq "T"){$Ts += $match;}
 
    print "$chr\t$pos\t$ref\t$covseqs\t$match\t$mismatch\t$As\t$Cs\t$Gs\t$Ts\n";
  }
}
close PILE;
