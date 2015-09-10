#!/usr/bin/perl

use strict;
use warnings;

my$const = $ARGV[0];
my$length = $ARGV[1];
my$output = $ARGV[2];

my$chr; my$pos; my$exon; my$start; my$stop; my$begin; my$end; my$max;
my@elements;
my%regs;

open(CONST,"$const") or die "\nError opening $const\n";
while(<CONST>)
{
  chomp;
  @elements = split("\t",$_);

  # mel_mel
  $chr = $elements[0];
  $start = $elements[1];
  $stop = $elements[2];
  $exon = $elements[3];
  if($stop-$start >= $length){$regs{$chr}{$exon} = [$start,$stop];}


  # mel_sim
  #$chr = $elements[0];
  #$pos = $elements[1]; # 1-based position in SNP BED file since 1-based in pileup
  #$base = $elements[3]; # if reference genome was used
  #$base = $elements[2]; # if alternative genome was used

}
close CONST;

open(OUT,"> $output") or die "Error writing to $output\n";

foreach $chr (keys %regs)
{
  foreach $exon (keys %{$regs{$chr}})
  {
    $begin = $regs{$chr}{$exon}[0];
    $max = $regs{$chr}{$exon}[1] - $length;
    
    while($begin+$length <= $max)
    {
      $end = $begin + $length;
      print OUT "$chr\t$begin\t$end\n";
      $begin += 1;
    }
  }
}
close OUT;
