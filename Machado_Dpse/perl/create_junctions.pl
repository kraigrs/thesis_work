#!/usr/bin/perl

######################################################################################
# 
# 03/21/2013
#
# create_transcripts.pl
#
# Purpose: using a list of exons, create every possible junction in a FASTA format
#
######################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$list = $ARGV[1];

my$chr; my$start; my$stop; my$name; my$coords; my$gene; my$line; my$strand; my$i; my$j;
my$e1; my$e2; my$e1_start; my$e2_start; my$e1_stop; my$e2_stop; my$string1; my$string2; my$length1; my$length2;
my$first = 1; my$seq = "";
my@elements; my@meta; my@regs;
my%data; my%genome;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome{$chr} = $seq;} # the first time, don't save sequence in hash
    $chr = $1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $chr\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

# read in exons
open(BED,"$list") or die "\nError opening $list\n";
while(<BED>)    
{
  chomp;
  $line = $_;
  @elements = split(/\s+/,$line);
  $chr = $elements[0];
  $start = $elements[1];
  $stop = $elements[2];
  $name = $elements[3];
  #$strand = $elements[4];

  @meta = split(/\_/,$name);
  $gene = $meta[0];

  #$coords = $start."_".$stop;
  #push @{ $data{$chr}{$gene} }, $coords;

  $data{$chr}{$gene}{$start} = $stop;
}
close BED;

foreach $chr (keys %data)
{
  foreach $gene (keys %{$data{$chr}})
  {
    if($genome{$chr})
    {
      foreach $start (sort {$a<=>$b} keys %{ $data{$chr}{$gene} })
      {
        push @regs, $start."_".$data{$chr}{$gene}{$start};
      }
      if(scalar(@regs) > 1)
      {
        for($i=0;$i<@regs;$i++)
        {
          for($j=0;$j<@regs;$j++)
	  {
            if($i != $j && $i < $j)
            {
              #print "i:$i\tj:$j\n";

              $e1 = $regs[$i];
              @meta = split(/\_/,$e1);
              $e1_start = $meta[0];
              $e1_stop = $meta[1];

              $e2 = $regs[$j];
              @meta = split(/\_/,$e2);
              $e2_start = $meta[0];
              $e2_stop = $meta[1];

              if($e1_stop <= $e2_start)
	      {
                $length1 = $e1_stop-$e1_start;
                $string1 = substr $genome{$chr}, $e1_start, $length1;

                $length2 = $e2_stop-$e2_start;
                $string2 = substr $genome{$chr}, $e2_start, $length2;

                print "\>$chr\_$gene\_$e1\:$e2\n";
                print $string1.$string2."\n";
	      }
            }
          }
	}
      }
    }
    @regs = ();
  }
}
