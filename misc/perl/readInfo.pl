#!/usr/bin/perl

###################################################################################
#
# 11/05/2009
#
# readInfo.pl
#
# Purpose: figure out how many sequences cover each exon in a gene (attempt at 
#          elucidating constitutive exons in genes) --> also try and find the 
#          sequences that  may contain more than 1 SNP (still trying to figure
#          that one out)
# 
# Input: a SAM alignment file from Bowtie will be used for the short reads.
#        First, a large hash containing all the genes, their exons, and start and
#        end coordinates for said exons will be made. Then we will parse through 
#        the SAM alignment file from Bowtie and, for each read, be able to tell
#        which exon that read falls in
#
# Output: a list of all the genes and their exons with counts of UNIQUE sequences
#         that are contained within each exon, and which species they align to 
#
###################################################################################

use strict;
use warnings;

my$starttime = time;

my$line;
my@elements;
my$gene;
my@loc;
my@exons;
my$i;
my$exon;
my$start;
my$stop;
my$string;
my@array;
my@AoA;
my%fastaHash;

my$gene_ref = $ARGV[0];
# open(GENES,"$gene_ref") or die "\nError opening $gene_ref\n";
# while(<GENES>)
# {
#   chomp $_;  
#   if(/^>(.+)\s*/)
#   { 
#     $line = $1;
#     @elements = split(/\s/,$line);
#     $gene = $elements[0]; 
#     #print "$gene\n";
#     #print "$elements[2]\n";
#     @loc = split(/\,/,$elements[2]);

#     foreach(@loc)
#     {
#       if(/(\d+\.\.\d+)/)
#       {
#         #print "$1\n"; 
#         push @exons,$1;
#       }
#     }
#     #print "@exons\n";

#     for($i = 0; $i < @exons; $i++)
#     {
#       $exon = $i + 1;
#       $fastaHash{$gene}{$exon} = $exons[$i];
#       #print "$gene\t$exon\t$fastaHash{$gene}{$exon}\n";
#     }
#     splice(@exons);
#   }
# }

# close GENES;
# print "\nFinished making gene hash!\n";


#print "gene\texon\tposition";
#for my$k1 (sort keys %fastaHash)
#{
#  for my$k2 (sort {$a <=> $b} keys %{$fastaHash{$k1}})
#  {
#    print "\n$k1\t$k2\t$fastaHash{$k1}{$k2}\n";
#  }
#}

my$pos;
my$ref;
my$cov;
my$call;
my@alleles;
my$match;
my$As;
my$Cs;
my$Gs;
my$Ts;
my$Ns;
my%SNPhash;
my%geneCov;

my$SNP_file = $ARGV[1];
open(SNP,"$SNP_file") or die "\nError opening $SNP_file\n";
while(<SNP>)
{
  $match = 0;
  $As = 0;
  $Cs = 0;
  $Gs = 0;
  $Ts = 0; 
  $Ns = 0;

  chomp $_;  
  $line = $_;
  @elements = split(/\s/,$line);
  $gene = $elements[0];
  $pos = $elements[1];
  $ref = $elements[2];
  $cov = $elements[7];
  $call = $elements[8];
  #print "Call:{$call}\t";
  @alleles = split("",$call);
  #foreach(@alleles){print "$_";}

  foreach(@alleles)
  {
    if(/[.,]/){$match+=1;}
    elsif(/[Aa]/){$As+=1;}
    elsif(/[Cc]/){$Cs+=1;}
    elsif(/[Gg]/){$Gs+=1;}
    elsif(/[Tt]/){$Ts+=1;}
    elsif(/[Nn]/){$Ns+=1;}
  }
  #print "$gene\t$pos\t$call\t$As\t$Cs\t$Gs\t$Ts\n";

  unless($cov<5 && (($As>0&&$Cs>0)||($As>0&&$Gs>0)||($As>0&&$Ts>0)||($Cs>0&&$Gs>0)||($Cs>0&&$Ts>0)||($Gs>0&&$Ts>0)||$Ns>0))
  {
    unless($match==$cov)
    {
      $SNPhash{$gene}{$pos} = $ref;
      #$geneCov{$gene} += $cov;
    }
  }    
}
close SNP;
print "\nFinished making SNP hash!\n";

#print "gene\tpos\tconsensus\n";
#for my$k1 (sort keys %SNPhash)
#{
#  for my$k2 (sort keys %{$SNPhash{$k1}})
#  {
#    print "$k1\t$k2\t$SNPhash{$k1}{$k2}\n";
#  }
#}

my$mate1found = 0;
my$mate2found = 0;
my$bothfound = 0;
my$read1;
my$read2;
my$gene1;
my$gene2;
my$m1start;
my$m2start;
my$seq1;
my$seq2;
my@seq1;
my@seq2;
my$exonKey;
my$exonNumber;
my$exonHits;
my$posKey;
my%geneHash;

my$mel1Ct;
my$sim1Ct;
my$mel2Ct;
my$sim2Ct;

my$m;
my$s;
my$b;
my$flag;
my$MEL;
my$SIM;

my$matepairs = $ARGV[2];
my$read_out = $ARGV[3];
print "\nParsing SAM output...\n";
open(SAM,"$matepairs") or die "\nError opening $matepairs\n";
while(<SAM>)
{
  chomp $_;  
  if(/^HWI.+\/1.+/)
  {
    $mate1found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read1 = $1;}
    #print "mate 1 found: $read1\n";
    $gene1 = $elements[2];
    $m1start = $elements[3]; 
    $seq1 = $elements[9];   
  }  
  elsif(/^HWI.+\/2.+/)
  {
    $mate2found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read2 = $1;}
    #print "mate 2 found: $read2\n";
    $gene2 = $elements[2];
    $m2start = $elements[3]; 
    $seq2 = $elements[9];
  }

  if($mate1found == 1 && $mate2found == 1 && $read1 eq $read2 && $gene1 eq $gene2)
  {
    #print "$read1\n";
    #print "Found both mates! Entering loop...\n";
    #print "$gene1\n";
    $mel1Ct = 0;
    $sim1Ct = 0;
    $mel2Ct = 0;
    $sim2Ct = 0;
    $m = 0;
    $s = 0;
    $b = 0;
    $flag = 'no-call';
    $MEL = 0;
    $SIM = 0;
    @seq1 = split("",$seq1);
    #print "@seq1\n";
    @seq2 = split("",$seq2);
    #print "@seq2\n";
    for $posKey (keys %{$SNPhash{$gene1}}) 
    {
      #print "$posKey\n";
      if($posKey >= $m1start && $posKey <= ($m1start+scalar@seq1-1))
      {
        if($seq1[$posKey-$m1start] eq $SNPhash{$gene1}{$posKey}){$mel1Ct += 1;}
        else{$sim1Ct += 1;}
      }
      elsif($posKey >= $m2start && $posKey <= ($m2start+scalar@seq2-1))
      {
        if($seq2[$posKey-$m2start] eq $SNPhash{$gene1}{$posKey}){$mel2Ct += 1;}
        else{$sim2Ct += 1;} 
      } 
    }

    $m = $mel1Ct + $mel2Ct;   
    $s = $sim1Ct + $sim2Ct;  

    if($m + $s == 0)
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;  
      $flag = 'no_SNPs';
      $b = 1;
    }
    elsif($m > 0 && $s == 0)
    {
      $geneHash{$gene1}{'m'} += 1; 
      $geneHash{$gene1}{'s'} += 0;
      $flag = 'confident';
      $MEL = 1;
    }
    elsif($m == 0 && $s > 0)
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 1;     
      $flag = 'confident';
      $SIM = 1;
    }
    elsif(($mel1Ct > 1 && $sim1Ct == 1) && ($mel2Ct == 0 && $sim2Ct == 0))
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;  
      $flag = 'sequencing_error';
    }
    elsif(($mel1Ct == 0 && $sim1Ct == 0) && ($mel2Ct > 1 && $sim2Ct == 1))
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;  
      $flag = 'sequencing_error';
    } 
    elsif(($mel1Ct == 1 && $sim1Ct > 1) && ($mel2Ct == 0 && $sim2Ct == 0))
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;  
      $flag = 'sequencing_error';
    }
    elsif(($mel1Ct == 0 && $sim1Ct == 0) && ($mel2Ct == 1 && $sim2Ct > 1))
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;  
      $flag = 'sequencing_error';
    }
    elsif(($mel1Ct >= 1 && $sim1Ct == 0) && ($mel2Ct == 0 && $sim2Ct >= 1))
    {
      $geneHash{$gene1}{'m'} += 1; 
      $geneHash{$gene1}{'s'} += 1; 
      $flag = 'trans-splicing';
      $MEL = 1;
      $SIM = 1;
    }
    elsif(($mel1Ct == 0  && $sim1Ct >= 1) && ($mel2Ct >= 1 && $sim2Ct == 0))
    {
      $geneHash{$gene1}{'m'} += 1;
      $geneHash{$gene1}{'s'} += 1; 
      $flag = 'trans-splicing';
      $MEL = 1;
      $SIM = 1;
    }
    elsif( (($mel1Ct > 0 && $sim1Ct < 2) && ($mel2Ct > 0 && $sim2Ct == 0)) || (($mel1Ct > 0 && $sim1Ct == 0) && ($mel2Ct > 0 && $sim2Ct < 2)) )
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;
      $flag = 'sequencing_error';
    }
    elsif( (($mel1Ct < 2 && $sim1Ct > 0) && ($mel2Ct == 0 && $sim2Ct > 0)) || (($mel1Ct == 0 && $sim1Ct > 0) && ($mel2Ct <= 1 && $sim2Ct > 0)) )
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;
      $flag = 'sequencing_error';
    }
    else
    {
      $geneHash{$gene1}{'m'} += 0;
      $geneHash{$gene1}{'s'} += 0;
      $flag = 'no_call';
    }
    $geneHash{$gene1}{'b'} += $b;
   
    open(OUT,">> $read_out") or die "Error writing to $read_out\n";
    print OUT "\n$read1\t$gene1\t$seq1\t$seq2\t$mel1Ct\t$sim1Ct\t$mel2Ct\t$sim2Ct\t$b\t$MEL\t$SIM\t$flag";
    close OUT;

    $mate1found = 0;
    $mate2found = 0;
  }
}
close SAM;

open(SAM,"$singlereads") or die "\nError opening $singlereads\n";
while(<SAM>)
{
  chomp $_;  
  if(/^HWI.+\/1.+/)
  {
    $mate1found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read1 = $1;}
    #print "mate 1 found: $read1\n";
    $gene1 = $elements[2];
    $m1start = $elements[3]; 
    $seq1 = $elements[9];   
  }  
  elsif(/^HWI.+\/2.+/)
  {
    $mate2found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read2 = $1;}
    #print "mate 2 found: $read2\n";
    $gene2 = $elements[2];
    $m2start = $elements[3]; 
    $seq2 = $elements[9];
  }  

  if($mate1found == 1)
  {
    $mel1Ct = 0;
    $sim1Ct = 0;
    $m = 0;
    $s = 0;
    $b = 0;
    $flag = 'no-call';
    $MEL = 0;
    $SIM = 0;
    @seq1 = split("",$seq1);
    #print "@seq1\n";
    for $posKey (keys %{$SNPhash{$gene1}}) 
    {
      if($posKey >= $m1start && $posKey <= ($m1start+scalar@seq1-1))
      {
        if($seq1[$posKey-$m1start] eq $SNPhash{$gene1}{$posKey}){$mel1Ct += 1;}
        else{$sim1Ct += 1;}
      }
    }
    if($mel1Ct + $sim1Ct == 0)
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 0;  
      $flag = 'no_SNPs';
      $b = 1;
    }
    elsif($mel1Ct > 0 && $sim1Ct == 0)
    {
      $geneHash{$gene1}{'m'} += 1; 
      $geneHash{$gene1}{'s'} += 0;
      $flag = 'confident';
      $MEL = 1;
    }
    elsif($mel1Ct == 0 && $sim1Ct > 0)
    {
      $geneHash{$gene1}{'m'} += 0; 
      $geneHash{$gene1}{'s'} += 1;     
      $flag = 'confident';
      $SIM = 1;
    }
    open(OUT,">> $read_out") or die "Error writing to $read_out\n";
    print OUT "\n$read1\t$gene1\t$seq1\tNA\t$mel1Ct\t$sim1Ct\tNA\tNA\t$b\t$MEL\t$SIM\t$flag";
    close OUT;
    
    $mate1found = 0;
  }
  elsif($mate2found == 1)
  {
    $mel2Ct = 0;
    $sim2Ct = 0;
    $m = 0;
    $s = 0;
    $b = 0;
    $flag = 'no-call';
    $MEL = 0;
    $SIM = 0;
    for $posKey (keys %{$SNPhash{$gene2}}) 
    {
      #print "$posKey\n";
      if($posKey >= $m2start && $posKey <= ($m2start+scalar@seq2-1))
      {
        if($seq2[$posKey-$m2start] eq $SNPhash{$gene2}{$posKey}){$mel2Ct += 1;}
        else{$sim2Ct += 1;} 
      } 
    }
    if($mel2Ct + $sim2Ct == 0)
    {
      $geneHash{$gene2}{'m'} += 0; 
      $geneHash{$gene2}{'s'} += 0;  
      $flag = 'no_SNPs';
      $b = 1;
    }
    elsif($mel2Ct > 0 && $sim2Ct == 0)
    {
      $geneHash{$gene2}{'m'} += 1; 
      $geneHash{$gene2}{'s'} += 0;
      $flag = 'confident';
      $MEL = 1;
    }
    elsif($mel2Ct == 0 && $sim2Ct > 0)
    {
      $geneHash{$gene2}{'m'} += 0; 
      $geneHash{$gene2}{'s'} += 1;     
      $flag = 'confident';
      $SIM = 1;
    }
    open(OUT,">> $read_out") or die "Error writing to $read_out\n";
    print OUT "\n$read2\t$gene2\tNA\t$seq2\tNA\tNA\t$mel2Ct\t$sim2Ct\t$b\t$MEL\t$SIM\t$flag";
    close OUT;

    $mate2found = 0;
  }
}
close SAM;

my$madj;
my$sadj;

my$gene_out = $ARGV[4];
open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
print OUT "gene\tmel\tsim\tboth\tmelAdj\tsimAdj";
close OUT;

for my$k1 (sort keys %geneHash)
{
  if(($geneHash{$k1}{'m'} > 0 || $geneHash{$k1}{'s'} > 0) && $k1 ne "*")
  {
    $madj = $geneHash{$k1}{'m'} + ($geneHash{$k1}{'m'} / ($geneHash{$k1}{'m'} + $geneHash{$k1}{'s'})) * $geneHash{$k1}{'b'};  
    $sadj = $geneHash{$k1}{'s'} + ($geneHash{$k1}{'s'} / ($geneHash{$k1}{'m'} + $geneHash{$k1}{'s'})) * $geneHash{$k1}{'b'};
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$k1\t$geneHash{$k1}{m}\t$geneHash{$k1}{s}\t$geneHash{$k1}{b}\t$madj\t$sadj";
    close OUT;
  }
}

print "Done!\n";
printf ("\nTime elapsed: %d\n\n",time-$starttime);

#         for $exonKey (sort {$a <=> $b} keys %{$fastaHash{$gene1}})
#         {
#           if($exonKey == 1)
#           {
#             $start = 1;
#             $stop = (split(/\.\./,$fastaHash{$gene1}{$exonKey}))[1] - (split(/\.\./,$fastaHash{$gene1}{$exonKey}))[0] + $start;  
          
#             if(($m1start >= $start) && ($m1start+scalar@seq1-1 <= $stop))
#             {
#               $exonNumber = $exonKey;
#               #print "$exonNumber\n";
#               $exonHits+=1;
#             } 
#           }
#           else
#           {
#             $start = $stop + 1;
#             $stop = (split(/\.\./,$fastaHash{$gene1}{$exonKey}))[1] - (split(/\.\./,$fastaHash{$gene1}{$exonKey}))[0] + $start;         
              
#             if(($m1start >= $start) && ($m1start+scalar@seq1-1 <= $stop))
#             {
#               $exonNumber = $exonKey;
#               #print "$exonNumber\n";
#               $exonHits+=1;
#             }  
#           }  
#         }




















