#!/usr/bin/perl

######################################################################################
# 
# 04/27/2012
#
# get_mappability_exons.pl
#
# Purpose: using BED file of exons, extract mappability measures from GEM output
# 
# Input: BED file of exons (or other regions I guess), GEM mappability output from FASTA genome        
#
# Output: for each input in BED file, output a mappability measure (average? all values? not sure yet...)
#
######################################################################################

use strict;
use warnings;

my$BED = $ARGV[0];
my$GEM = $ARGV[1];
my$output = $ARGV[2];

my$encode = 0; my$fasta = 0; my$counter = 0; my$seq = ""; my$first = 1; my$readlen = 0;
our($chr,$start,$stop,$gene,$locus,$code,$freq,$begin,$end,$len,$string,$aref,$mappability,$kmer);
our(@elements,@array);
our(%regs,%map,%genome);

open(BED,"$BED") or die "\nError opening $BED\n";
while(<BED>)    
{
  chomp;
  @elements = split(/\s+/,$_);
  $chr = $elements[0]; 
  $start = $elements[1];
  $stop = $elements[2];
  $gene = $elements[3];
  $locus = $gene."_".$start."_".$stop;

  $regs{$chr}{$locus} = [$start,$stop];
  #$regs{$chr}{$gene} = [$start,$stop];
}
close BED;

open(GEM,"$GEM") or die "\nError opening $GEM\n";
while(<GEM>)    
{
  chomp;
  #$counter += 1;
  #print "$counter\n";

  if(/READ LENGTH/)
  {
    $readlen = 1;
  }
  elsif($readlen == 1)
  {
    if(/^(\d+)$/)
    {
      $kmer = $1;
      $readlen = 0;
    }
  }

  if(/ENCODING/)
  {
    $encode = 1;
  }
  elsif($encode == 1)
  {
    if(/\'(\S{1})\'\~\[(\d+)\-\d+\]/)
    {
      $code = $1; $freq = $2;
      #print "\nCode: $code\tFreq: $freq\n";
      $map{$code} = $freq;
    }
  }

  if(/^\~(chr\S+)$/)
  {
    unless($first == 1){$genome{$chr} = $seq;}

    $chr = $1;
    $seq = "";
    #print "\nChromosome: $chr\n";

    $fasta = 1;
    $encode = 0;
    $first = 0;
  }
  elsif($fasta == 1)
  {
    $seq = $seq.$_;
    #print "\nSequence: $seq\n";
  }
}
close GEM;
$genome{$chr} = $seq;

#print "\nDone\n"; exit;

foreach $chr (keys %genome)
{
  $len = length($genome{$chr});
  print "\nChromosome: $chr\tLength: $len\n";
}
print "\n";

open(OUT,">$output") or die "Error writing to $output\n";

#print OUT "locus\tsum\tlength\tmappability\n";
print OUT "locus\tsum\tlength";

foreach $chr (keys %regs)
{
  foreach $locus (keys %{$regs{$chr}})
  {
    if($genome{$chr})
    {
      $begin = $regs{$chr}{$locus}[0];
      #$end = $regs{$chr}{$locus}[1]-$regs{$chr}{$locus}[0];
      #print "Kmer: $kmer\n";
      $end = $regs{$chr}{$locus}[1]-$regs{$chr}{$locus}[0]-$kmer+1; # mappability measured to last 
                                                                    # possible position in exon of read length k
      #print "Chromosome: $chr\n$genome{$chr}"; exit;

      unless($end < 1) # forgot initially to account for exons that are shorter than read length
      {
        $string = substr $genome{$chr}, $begin, $end;
        #print "Begin: $begin\tLength: $end\tString: $string";
        @array = split("",$string);
        $aref = \@array;
        $mappability = &aggregate($aref); # this subroutine simply sums

        #print "\nChromosome: $chr\tLocus: $locus\tMappability: $mappability\n";
        #print "@array\n";

        #print OUT "$locus\t$mappability\t$end\t$string\n";
        print OUT "\n$locus\t$mappability\t$end";
      }
    }
  }
}
close OUT;

sub aggregate
{
  my$ref = shift;
  my@a = @{$ref};
  my$sum = 0;
  my$ct = 0;  

  foreach(@a)
  {
    unless(!$map{$_})
    {
      $ct += 1;
      #print "Code: $_\tFreq: $map{$_}\n";
      #$sum += $map{$_};
      $sum += 1/$map{$_};
    }
  }
  return $sum;
}
