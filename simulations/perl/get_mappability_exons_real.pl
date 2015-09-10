#!/usr/bin/perl

######################################################################################
# 
# 01/17/2012
#
# get_mappability_exons_real.pl
#
######################################################################################

use strict;
use warnings;

my$GEM = $ARGV[0];
my$output = $ARGV[1];

my$encode = 0; my$fasta = 0; my$counter = 0; my$seq = ""; my$first = 1; my$readlen = 0; my$i; my$j;
our($chr,$start,$stop,$gene,$locus,$code,$freq,$begin,$end,$len,$string,$aref,$mappability,$kmer);
our(@elements,@array);
our(%map,%genome);

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
    #print "Here\n";

    if(/\'(.+)\'\~\[(\d+)\-\d+\]/)
    {
      $code = $1; $freq = $2;
      #print "\nCode: $code\tFreq: $freq\n";
      $map{$code} = $freq;
    }
  }

  if(/^\~{1}([^\~]+)$/)
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

#foreach $code (keys %map)
#{
#  print "\n\$map\{$code\} -> $map{$code}";
#}
#print "\n\n"; exit;



open(OUT,">$output") or die "Error writing to $output\n";

foreach $chr (keys %genome)
{
  $string = $genome{$chr};
  @array = split("",$string);
  for($i=0;$i<@array;$i++)
  {
    unless($array[$i] eq " ")
    {
      if($map{$array[$i]} != 1)
      {
        $j = $i+$kmer;
        print OUT "$chr\t$i\t$j\n";
      }
    }
  }
}
close OUT;
