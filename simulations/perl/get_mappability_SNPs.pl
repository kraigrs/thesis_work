#!/usr/bin/perl

######################################################################################
# 
# 04/27/2012
#
# get_mappability.pl
#
# Purpose: using BED file of SNPs, extract mappability measures from GEM output
# 
# Input: BED file of SNPs, GEM mappability output from FASTA genome        
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
our($chr,$start,$stop,$gene,$locus,$code,$freq,$begin,$end,$len,$string,$aref,$mappability,$kmer,$SNP,$temp);
our(@elements,@meta,@array,@return);
our(%regs,%map,%genome);

open(BED,"$BED") or die "\nError opening $BED\n";
while(<BED>)    
{
  chomp;
  @elements = split(/\s+/,$_);
  $chr = $elements[0]; 
  $SNP = $elements[2]; # 1-based SNP coordinate
  #$SNP = $elements[1]; # 1-based SNP coordinate

  #$start = $elements[5];
  #$stop = $elements[6];
  #$gene = $elements[7];
  #$locus = $gene."_".$start."_".$stop;

  #$regs{$chr}{$SNP} = [$start,$stop,$locus];
  $regs{$chr}{$SNP} = 1;
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
    #if(/\'(\S{1})\'\~\[(\d+)\-\d+\]/)
    if(/\'(.+)\'\~\[(\d+)\-\d+\]/)
    {
      $code = $1; $freq = $2;
      #print "\nCode: $code\tFreq: $freq\n";
      $map{$code} = $freq;
    }
  }

  if(/^\~(\S+)$/)
  {    
    $temp = $1;
    unless($first == 1){$genome{$chr} = $seq;}
    
    if($temp =~ /[\|]/)
    {
      @meta = split(/\|/,$temp);
      $chr = $meta[0];
    }
    else{$chr = $temp;}
    
    $seq = "";
    #print "\nChromosome: $chr\n";

    $fasta = 1;
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

#print "\nCode ' ': {$map{' '}}\n\n"; exit;
#print "\nDone\n"; exit;

#foreach $code (keys %map)
#{
#  print "\nCode '$code': {$map{$code}}\n";
#}
#exit;

#foreach $chr (keys %genome)
#{
#  $len = length($genome{$chr});
#  print "\nChromosome: $chr\tLength: $len\n";
#}
#exit;
#print "\n";

open(OUT,">$output") or die "Error writing to $output\n";

#print OUT "chr\tposition\tsum\tlength\tlocus";
print OUT "chr\tposition\tsum\tlength";

foreach $chr (keys %regs)
{
  foreach $SNP (keys %{$regs{$chr}})
  {
    #print "Chrom: $chr\tSNP: $SNP\n";
    #print "\{$genome{$chr}\}"; exit;

    if($genome{$chr})
    {
      #if($chr eq "F10014_SI"){print "$chr\n$genome{$chr}"; exit;}

      #$start = $regs{$chr}{$SNP}[0];
      #$stop = $regs{$chr}{$SNP}[1];
      #$locus = $regs{$chr}{$SNP}[2];

      $start = 0;
      $stop = length($genome{$chr});

      unless($stop - $start < $kmer)
      {
        if($SNP - $kmer < $start){$begin = $start;} # if SNP is too close to beginning of exon
        else{$begin = $SNP - $kmer;}

        if($SNP + $kmer - 1 > $stop){$end = $stop - $SNP + 1;} # if SNP is too close to end of exon
        else{$end = $SNP - $begin;}

        #print "$chr\t$SNP\t$start\t$stop\t$begin\t$end\n";

        $string = substr $genome{$chr}, $begin, $end;
        @array = split("",$string);
        $aref = \@array;
        @return = &aggregate($aref); # this subroutine simply sums up the mappabilities

        #print OUT "\n$chr\t$SNP\t$return[0]\t$return[1]\t$locus";
        print OUT "\n$chr\t$SNP\t$return[0]\t$return[1]";
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
    #unless(!$map{$_})
    unless($map{$_} == 0)
    {
      $ct += 1;
      #print "Code: $_\tFreq: $map{$_}\n";
      #$sum += $map{$_};
      $sum += 1/$map{$_};
    }
  }
  return ($sum,$ct);
}
