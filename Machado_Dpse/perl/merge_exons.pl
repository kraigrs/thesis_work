#!/usr/bin/perl

# for genes with multiple exons, merge them using mergeBed

use strict;
use warnings;

my$file = $ARGV[0];
my$out = $ARGV[1];

my$line; my$chr; my$start; my$stop; my$gene; my$loc; my$strand; my$exon; my$entry; my$i;
my@elements; my@meta; my@coords;
my%geneHash; my%orient;

open(BED,"$file") or die "\nError opening $file\n";
while(<BED>)    
{
  chomp;
  $line = $_;

  if(/^track/){next;}
  else
  {
    @elements = split("\t",$line);

    $chr = $elements[0]; 
    $start = $elements[1];
    $stop = $elements[2];
    $loc = $elements[3];
    $strand = $elements[5];

    @meta = split(/\_/,$loc);
    $gene = $meta[0];
    push @{ $geneHash{$gene} }, "$line";
    $orient{$gene} = $strand;
    #print "$line\n$geneHash{$gene}[0]\n"; exit;
  }
}
close BED;

#$count = keys %geneHash;
#print "\n$count\n\n";

open(OUT,"> $out") or die "\nError writing to $out\n";

foreach $gene (keys %geneHash)
{
  if(scalar(@{ $geneHash{$gene} }) == 1)
  {
    @elements = split("\t",$geneHash{$gene}[0]);

    $chr = $elements[0]; 
    $start = $elements[1];
    $stop = $elements[2];
    $loc = $elements[3];
    $strand = $elements[5];

    @meta = split(/\_/,$loc);
    print OUT "$chr\t$start\t$stop\t$meta[0]\:$meta[1]\_$meta[2]\t0\t$strand\n";
  }
  else
  {
    # for each gene with multiple exons, put those in a temporary file

    open(TEMP,"> temp") or die "\nError writing to temp\n";
    foreach $entry (@{ $geneHash{$gene} })
    {
      print TEMP "$entry\n";
    }
    close TEMP;
    
    system "mergeBed -nms -i temp > merged";
    #exit;

    open(MERGED,"merged") or die "\nError opening merged\n";
    while(<MERGED>)    
    {
      chomp;
      $line = $_;
      @elements = split("\t",$line);
      $chr = $elements[0];
      $start = $elements[1];
      $stop = $elements[2];
      $loc = $elements[3];

      if($loc =~ /^(CG\d+)\_\S+/){$gene = $1;}
      #print "$line\n$gene\n";
  
      @meta = split(/\;/,$loc);
      foreach $exon (@meta)
      {
        if($exon =~ /^\S+\_(\d+\_\d+)$/){push @coords, $1;}
      }
      print OUT "$chr\t$start\t$stop\t$gene:";
      for($i=0;$i<@coords;$i++)
      {
        if($i == scalar(@coords)-1){print OUT "$coords[$i]";}
        else{print OUT "$coords[$i],";}
      }
      print OUT "\t0\t$orient{$gene}\n";
      @coords = ();
    }
    close MERGED;
    #exit;
  }
}

