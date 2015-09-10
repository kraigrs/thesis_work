#!/usr/bin/perl

use strict;
use warnings;

my$bed = $ARGV[0];
#my$out_dir = $ARGV[1];

our($line,$chrA,$startA,$stopA,$geneA,$locA,$strandA,$exonA,$chrB,$startB,$stopB,$geneB,$locB,$strandB,$exonB,$key);
our(@elements,@stuff);
our(%same,%check);

system "intersectBed -a $bed -b $bed -wa -wb > $bed\.intersected";

#open(GENE,">$out_dir\/$bed\.gene") or die "\nError writing to $out_dir\/$bed\.gene\n";
#open(SAME,">$out_dir\/$bed\.same") or die "\nError writing to $out_dir\/$bed\.same\n";
#open(DIFF,">$out_dir\/$bed\.diff") or die "\nError writing to $out_dir\/$bed\.diff\n";
#open(NORM,">$out_dir\/$bed\.norm") or die "\nError writing to $out_dir\/$bed\.norm\n";

open(BED,"$bed\.intersected") or die "\nError opening $bed\.intersected\n";
while(<BED>)    
{
  chomp;
  $line = $_;

  if(/^track/){next;}
  else
  {
    @elements = split("\t",$line);

    $chrA = $elements[0]; 
    $startA = $elements[1];
    $stopA = $elements[2];
    $locA = $elements[3];
    #if($locA =~ /(\S+)\-R[A-Z]/){$geneA = $1;}
    @stuff = split(/\_/,$locA);
    $geneA = $stuff[0];
    $exonA = $stuff[1]."_".$stuff[2];
    $strandA = $elements[5];

    $chrB = $elements[6]; 
    $startB = $elements[7];
    $stopB = $elements[8];
    $locB = $elements[9];
    #if($locB =~ /(\S+)\-R[A-Z]/){$geneB = $1;}
    @stuff = split(/\_/,$locB);
    $geneB = $stuff[0];
    $exonB = $stuff[1]."_".$stuff[2];
    $strandB = $elements[11];

    if($chrA eq $chrB && $geneA ne $geneB && $strandA eq $strandB && $exonA eq $exonB)
    {
      # same coords, different gene names

      #print SAME "$line\n";
      print "$chrA\t$startA\t$stopA\t$locA\;$locB\n";
      #unless($same{$locA} && $same{$locB}){print SAME "$chrA\t$startA\t$stopA\t$locA\;$locB\t0\t$strandA\n";}
      #unless($same{$locA})
      #{
        #print SAME "$line\n"; 
        #print DIFF "$chrB\t$startB\t$stopB\t$locB\t0\t$strandB\t$chrB\t$startB\t$stopB\t$locB\t0\t$strandB\n"; 
        #$same{$locB} = $locA;
      #}
      #$check{$locA} = 1; $check{$locB} = 1;
    }

    #if($chrA eq $chrB && $geneA eq $geneB && $exonA ne $exonB && $strandA eq $strandB)
    #{
    #  # overlapping exons within a gene

    #  #print GENE "$line\n";
    #  print GENE "$chrA\t$startA\t$stopA\t$locA\t0\t$strandA\n";
    #  #print GENE "$chrA\t$startA\t$stopA\t$locA\t0\t$strandA\n";
    #  #$check{$locA} = 1; $check{$locB} = 1;
    #}

    #elsif($chrA eq $chrB && ($startA != $startB || $stopA != $stopB || $locA ne $locB))
    #{
    #  # overlaps between genes, not strand specific 

    #  print DIFF "$line\n";
    #  #unless($same{$locA} || $same{$locB}){print DIFF "$line\n";}
    #  #unless($check{$locA}){print DIFF "$chrA\t$startA\t$stopA\t$locA\t0\t$strandA\n";}
    #  #$check{$locA} = 1; $check{$locB} = 1;
    #}

    #elsif($chrA eq $chrB && $startA == $startB && $stopA == $stopB && $locA eq $locB && $strandA eq $strandB)
    #{
    #  # this should simply be the exons that don't have any overlapping info
    #  #unless($check{$locA} || $same{$locA}){print NORM "$chrA\t$startA\t$stopA\t$locA\t0\t$strandA\n";}
    #  print NORM "$line\n";
    #}
    #else
    #{
    #  # everything else
    #  print "$line\n";
    #}
  }
}
close BED;

#close DIFF; 
#close SAME; 
#close GENE; 
#close NORM;

#foreach $key (keys %overlap){print "$key\n";}
