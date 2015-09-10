#!/usr/bin/perl

########################################################################
# 
# 05/11/2010
#
# simIlluminaPE0.2.pl
#
# Purpose: simulate sets of paired-end reads from reference databases
# 
# Input: a FASTA reference and SNP locations (with bp info)
#
# Output: millions of short sequence reads with various properties (mutation rates,
#         lengths, etc.)
# 
# Fixed: new version 0.2 --> generate reads around SNPs (categorized as "both"), 
#                            Basically, generate reads from the entire reference
#
#        05/14/2010 --> created gene (or chrom) hash structure instead of going through file
#
# Notes: 05/13/2010 --> need to fix gene (or chromosome) loop so that we only
#                       go through the reference file once
#
########################################################################

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

srand(12345);

my$length = $ARGV[0]; # read length
my$insert = $ARGV[1]; # length between paired ends
my$error = $ARGV[2]; # base-calling error rate

my$match; my$As; my$Cs; my$Gs; my$Ts; my$Ns;
my$gene; my$pos; my$cons; my$cov; my$call;
my@elements; my@alleles;
my%SNPhash;

# now create a hash of SNPs to be used later to generate reads around
my$SNP_file = $ARGV[3];
open(SNP,"$SNP_file") or die "\nError opening $SNP_file\n";
while(<SNP>)
{
  $match = 0;
  $As = 0;
  $Cs = 0;
  $Gs = 0;
  $Ts = 0; 
  $Ns = 0;

  chomp;  
  @elements = split(/\s+/,$_);
  $gene = $elements[0];
  $pos = $elements[1];
  $cons = $elements[3];
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

  unless($cov<20 && (($As>0&&$Cs>0)||($As>0&&$Gs>0)||($As>0&&$Ts>0)||($Cs>0&&$Gs>0)||($Cs>0&&$Ts>0)||($Gs>0&&$Ts>0)||$Ns>0))
  {
    unless($match==$cov)
    {
      $SNPhash{$gene}{$pos} = $cons;
    }
  }    
}
close SNP;

my$seq; 
my%geneHash;

my$ref = $ARGV[4];
open(REF,"$ref") or die "\nError opening $ref\n";

while(<REF>)    
{
  chomp;
  if(/^>(\S+)/)
  {
    $gene = $1;
    undef $seq;
  }
  else
  {
    $seq = $seq . $_;
    $geneHash{$gene} = $seq;  
  }
}
close REF;

#foreach $gene (keys %geneHash)
#{
#  print "Gene: $gene\nSeq: $geneHash{$gene}\n";
#}

my$i; my$j; my$k = 1;
my$reflength; my$lengthTot; my$posNeg; my$false1; my$false2;
my$readRef1; my$readRef2; my$readCons1; my$readCons2;  
my$baseRef1; my$baseRef2; my$baseCons1; my$baseCons2;
my@nuc; my@nucNeg;

my$out1 = $ARGV[5];
my$out2 = $ARGV[6];

open(OUT1,">$out1") or die "\nError opening $out1\n";
open(OUT2,">$out2") or die "\nError opening $out2\n";

foreach $gene (keys %geneHash)
{
  @nuc = split("",$geneHash{$gene});
  @nucNeg = reverse @nuc;
  $reflength = scalar@nuc;
  $lengthTot = 2*$length+$insert; # total length of paired-end read
  #print "Length of reference: $reflength\n";

  foreach $pos (keys %{$SNPhash{$gene}})
  {
    unless($lengthTot > $reflength)
    {
      # positive orientation

      for($i=0;$i<$reflength-$lengthTot+1;$i++) # generate reads across entire reference
      {
        for($j=0;$j<$length;$j++)
        {
          if($error == 0) # no error model
	  {
            $readRef1 = $readRef1 . $nuc[$i+$j];
            $readRef2 = $readRef2 . $nuc[$i+$j+$length+$insert];
	    if($i+$j == $pos-1)
            {
              $readCons1 = $readCons1 . $SNPhash{$gene}{$pos};
              $readCons2 = $readCons2 . $nuc[$i+$j+$length+$insert];
            }
            elsif($i+$j+$length+$insert == $pos-1)
            {
              $readCons1 = $readCons1 . $nuc[$i+$j];
              $readCons2 = $readCons2 . $SNPhash{$gene}{$pos};
            }            
            else
            {
              $readCons1 = $readCons1 . $nuc[$i+$j];
              $readCons2 = $readCons2 . $nuc[$i+$j+$length+$insert];
            }
          }
          else # error model (much longer code!)
	  {
            $false1 = rand();
            $false2 = rand();

            #print "i = $i\tj = $j\n";
            if($false1 < $error && $false2 >= $error) # only error on mate 1
            {
              $baseRef1 = rand();
              $baseCons1 = rand();
 
              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef1: $baseRef1\tbaseCons1 = $baseCons1\n";
              if($baseRef1 < 0.25){$readRef1 = $readRef1 . "A";}
              elsif($baseRef1 >= 0.25 && $baseRef1 < 0.5){$readRef1 = $readRef1 . "C";}
              elsif($baseRef1 >= 0.5 && $baseRef1 < 0.75){$readRef1 = $readRef1 . "G";}
              else{$readRef1 = $readRef1 . "T";}

              if($baseCons1 < 0.25){$readCons1 = $readCons1 . "A";}
              elsif($baseCons1 >= 0.25 && $baseCons1 < 0.5){$readCons1 = $readCons1 . "C";}
              elsif($baseCons1 >= 0.5 && $baseCons1 < 0.75){$readCons1 = $readCons1 . "G";}
              else{$readCons1 = $readCons1 . "T";}
              
              $readRef2 = $readRef2 . $nuc[$i+$j+$length+$insert];
              if($i+$j+$length+$insert == $pos-1){$readCons2 = $readCons2 . $SNPhash{$gene}{$pos};}
              else{$readCons2 = $readCons2 . $nuc[$i+$j+$length+$insert];}
	    }
            elsif($false1 >= $error && $false2 < $error) # only error on mate 2
            {    
              $readRef1 = $readRef1 . $nuc[$i+$j+$length+$insert];
              if($i+$j+$length+$insert == $pos-1){$readCons1 = $readCons1 . $SNPhash{$gene}{$pos};}
              else{$readCons1 = $readCons1 . $nuc[$i+$j+$length+$insert];}

              $baseRef2 = rand();
              $baseCons2 = rand();

              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef2: $baseRef2\tbaseCons2 = $baseCons2\n";
              if($baseRef2 < 0.25){$readRef2 = $readRef2 . "A";}
              elsif($baseRef2 >= 0.25 && $baseRef2 < 0.5){$readRef2 = $readRef2 . "C";}
              elsif($baseRef2 >= 0.5 && $baseRef2 < 0.75){$readRef2 = $readRef2 . "G";}
              else{$readRef2 = $readRef2 . "T";}

              if($baseCons2 < 0.25){$readCons2 = $readCons2 . "A";}
              elsif($baseCons2 >= 0.25 && $baseCons2 < 0.5){$readCons2 = $readCons2 . "C";}
              elsif($baseCons2 >= 0.5 && $baseCons2 < 0.75){$readCons2 = $readCons2 . "G";}
              else{$readCons2 = $readCons2 . "T";}
	    }
            elsif($false1 >= $error && $false2 < $error) # error on mate 1 and mate 2
            {
              $baseRef1 = rand();              
              $baseRef2 = rand();
              $baseCons1 = rand();
              $baseCons2 = rand();
 
              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef1: $baseRef1\tbaseCons1 = $baseCons1\n";
              if($baseRef1 < 0.25){$readRef1 = $readRef1 . "A";}
              elsif($baseRef1 >= 0.25 && $baseRef1 < 0.5){$readRef1 = $readRef1 . "C";}
              elsif($baseRef1 >= 0.5 && $baseRef1 < 0.75){$readRef1 = $readRef1 . "G";}
              else{$readRef1 = $readRef1 . "T";}

              if($baseCons1 < 0.25){$readCons1 = $readCons1 . "A";}
              elsif($baseCons1 >= 0.25 && $baseCons1 < 0.5){$readCons1 = $readCons1 . "C";}
              elsif($baseCons1 >= 0.5 && $baseCons1 < 0.75){$readCons1 = $readCons1 . "G";}
              else{$readCons1 = $readCons1 . "T";}

              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef2: $baseRef2\tbaseCons2 = $baseCons2\n";
              if($baseRef2 < 0.25){$readRef2 = $readRef2 . "A";}
              elsif($baseRef2 >= 0.25 && $baseRef2 < 0.5){$readRef2 = $readRef2 . "C";}
              elsif($baseRef2 >= 0.5 && $baseRef2 < 0.75){$readRef2 = $readRef2 . "G";}
              else{$readRef2 = $readRef2 . "T";}

              if($baseCons2 < 0.25){$readCons2 = $readCons2 . "A";}
              elsif($baseCons2 >= 0.25 && $baseCons2 < 0.5){$readCons2 = $readCons2 . "C";}
              elsif($baseCons2 >= 0.5 && $baseCons2 < 0.75){$readCons2 = $readCons2 . "G";}
              else{$readCons2 = $readCons2 . "T";}
	    }
            else
	    {
              $readRef1 = $readRef1 . $nuc[$i+$j];
              $readRef2 = $readRef2 . $nuc[$i+$j+$length+$insert];
	      if($i+$j == $pos-1)
              {
                $readCons1 = $readCons1 . $SNPhash{$gene}{$pos};
                $readCons2 = $readCons2 . $nuc[$i+$j+$length+$insert];
              }
              elsif($i+$j+$length+$insert == $pos-1)
              {
                $readCons1 = $readCons1 . $nuc[$i+$j];
                $readCons2 = $readCons2 . $SNPhash{$gene}{$pos};
              }            
              else
              {
                $readCons1 = $readCons1 . $nuc[$i+$j];
                $readCons2 = $readCons2 . $nuc[$i+$j+$length+$insert];
              }
            }
	  }
	}
        if( ($pos-1 >= $i && $pos-1 <= $i+$length-1) || ($pos-1 >= $i+$length+$insert && $pos-1 <= $i+$lengthTot-1) ) # these are reads containing SNPs
	{
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\+\)\/1\n$readRef1\n";
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\+\)\/1\n$readCons1\n";

          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\+\)\/2\n$readRef2\n";
          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\+\)\/2\n$readCons2\n";
        }
        else # these are reads NOT containing SNPs
	{
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\+\)\/1\n$readRef1\n";
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\+\)\/1\n$readCons1\n";

          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\+\)\/2\n$readRef2\n";
          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\+\)\/2\n$readCons2\n";
        }

        undef $readRef1;
        undef $readRef2;
        undef $readCons1;
        undef $readCons2;

        $k += 1;    
      }
      $k = 1;
     
      # negative orientation

      $posNeg = $reflength - $pos + 1; # reverse the position of the SNP
      for($i=0;$i<$reflength-$lengthTot+1;$i++) # account for read exceeding reference
      {
        for($j=0;$j<$length;$j++)
        {
          if($error == 0) # no error model
	  {
            $readRef1 = $readRef1 . $nucNeg[$i+$j];
            $readRef2 = $readRef2 . $nucNeg[$i+$j+$length+$insert];
	    if($i+$j == $posNeg-1)
            {
              $readCons1 = $readCons1 . $SNPhash{$gene}{$pos};
              $readCons2 = $readCons2 . $nucNeg[$i+$j+$length+$insert];
            }
            elsif($i+$j+$length+$insert == $posNeg-1)
            {
              $readCons1 = $readCons1 . $nucNeg[$i+$j];
              $readCons2 = $readCons2 . $SNPhash{$gene}{$pos};
            }            
            else
            {
              $readCons1 = $readCons1 . $nucNeg[$i+$j];
              $readCons2 = $readCons2 . $nucNeg[$i+$j+$length+$insert];
            }
          }
          else # error model (much longer code!)
	  {
            $false1 = rand();
            $false2 = rand();

            #print "i = $i\tj = $j\n";
            if($false1 < $error && $false2 >= $error) # only error on mate 1
            {
              $baseRef1 = rand();
              $baseCons1 = rand();
 
              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef1: $baseRef1\tbaseCons1 = $baseCons1\n";
              if($baseRef1 < 0.25){$readRef1 = $readRef1 . "A";}
              elsif($baseRef1 >= 0.25 && $baseRef1 < 0.5){$readRef1 = $readRef1 . "C";}
              elsif($baseRef1 >= 0.5 && $baseRef1 < 0.75){$readRef1 = $readRef1 . "G";}
              else{$readRef1 = $readRef1 . "T";}

              if($baseCons1 < 0.25){$readCons1 = $readCons1 . "A";}
              elsif($baseCons1 >= 0.25 && $baseCons1 < 0.5){$readCons1 = $readCons1 . "C";}
              elsif($baseCons1 >= 0.5 && $baseCons1 < 0.75){$readCons1 = $readCons1 . "G";}
              else{$readCons1 = $readCons1 . "T";}
              
              $readRef2 = $readRef2 . $nucNeg[$i+$j+$length+$insert];
              if($i+$j+$length+$insert == $posNeg-1){$readCons2 = $readCons2 . $SNPhash{$gene}{$pos};}
              else{$readCons2 = $readCons2 . $nucNeg[$i+$j+$length+$insert];}
	    }
            elsif($false1 >= $error && $false2 < $error) # only error on mate 2
            {
              $readRef1 = $readRef1 . $nucNeg[$i+$j];
              if($i+$j == $posNeg-1){$readCons1 = $readCons1 . $SNPhash{$gene}{$pos};}
              else{$readCons1 = $readCons1 . $nucNeg[$i+$j];}

              $baseRef2 = rand();
              $baseCons2 = rand();

              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef2: $baseRef2\tbaseCons2 = $baseCons2\n";
              if($baseRef2 < 0.25){$readRef2 = $readRef2 . "A";}
              elsif($baseRef2 >= 0.25 && $baseRef2 < 0.5){$readRef2 = $readRef2 . "C";}
              elsif($baseRef2 >= 0.5 && $baseRef2 < 0.75){$readRef2 = $readRef2 . "G";}
              else{$readRef2 = $readRef2 . "T";}

              if($baseCons2 < 0.25){$readCons2 = $readCons2 . "A";}
              elsif($baseCons2 >= 0.25 && $baseCons2 < 0.5){$readCons2 = $readCons2 . "C";}
              elsif($baseCons2 >= 0.5 && $baseCons2 < 0.75){$readCons2 = $readCons2 . "G";}
              else{$readCons2 = $readCons2 . "T";}
	    }
            elsif($false1 < $error && $false2 < $error) # error on mate 1 and mate 2
            {
              $baseRef1 = rand();              
              $baseRef2 = rand();
              $baseCons1 = rand();
              $baseCons2 = rand();
 
              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef1: $baseRef1\tbaseCons1 = $baseCons1\n";
              if($baseRef1 < 0.25){$readRef1 = $readRef1 . "A";}
              elsif($baseRef1 >= 0.25 && $baseRef1 < 0.5){$readRef1 = $readRef1 . "C";}
              elsif($baseRef1 >= 0.5 && $baseRef1 < 0.75){$readRef1 = $readRef1 . "G";}
              else{$readRef1 = $readRef1 . "T";}

              if($baseCons1 < 0.25){$readCons1 = $readCons1 . "A";}
              elsif($baseCons1 >= 0.25 && $baseCons1 < 0.5){$readCons1 = $readCons1 . "C";}
              elsif($baseCons1 >= 0.5 && $baseCons1 < 0.75){$readCons1 = $readCons1 . "G";}
              else{$readCons1 = $readCons1 . "T";}

              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef2: $baseRef2\tbaseCons2 = $baseCons2\n";
              if($baseRef2 < 0.25){$readRef2 = $readRef2 . "A";}
              elsif($baseRef2 >= 0.25 && $baseRef2 < 0.5){$readRef2 = $readRef2 . "C";}
              elsif($baseRef2 >= 0.5 && $baseRef2 < 0.75){$readRef2 = $readRef2 . "G";}
              else{$readRef2 = $readRef2 . "T";}

              if($baseCons2 < 0.25){$readCons2 = $readCons2 . "A";}
              elsif($baseCons2 >= 0.25 && $baseCons2 < 0.5){$readCons2 = $readCons2 . "C";}
              elsif($baseCons2 >= 0.5 && $baseCons2 < 0.75){$readCons2 = $readCons2 . "G";}
              else{$readCons2 = $readCons2 . "T";}
	    }
            else
	    {
              $readRef1 = $readRef1 . $nucNeg[$i+$j];
              $readRef2 = $readRef2 . $nucNeg[$i+$j+$length+$insert];
	      if($i+$j == $posNeg-1)
              {
                $readCons1 = $readCons1 . $SNPhash{$gene}{$pos};
                $readCons2 = $readCons2 . $nucNeg[$i+$j+$length+$insert];
              }
              elsif($i+$j+$length+$insert == $posNeg-1)
              {
                $readCons1 = $readCons1 . $nucNeg[$i+$j];
                $readCons2 = $readCons2 . $SNPhash{$gene}{$pos};
              }            
              else
              {
                $readCons1 = $readCons1 . $nucNeg[$i+$j];
                $readCons2 = $readCons2 . $nucNeg[$i+$j+$length+$insert];
              }
            }
	  }
	}
        if( ($posNeg-1 >= $i && $posNeg-1 <= $i+$length-1) || ($posNeg-1 >= $i+$length+$insert && $posNeg-1 <= $i+$lengthTot-1) ) # these are reads containing SNPs
	{
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\-\)\/1\n$readRef1\n";
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\-\)\/1\n$readCons1\n";

          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\-\)\/2\n$readRef2\n";
          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\-\)\/2\n$readCons2\n";
        }
        else # these are reads NOT containing SNPs
	{
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\-\)\/1\n$readRef1\n";
          print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\-\)\/1\n$readCons1\n";

          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\-\)\/2\n$readRef2\n";
          print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\-\)\/2\n$readCons2\n";
        }

        undef $readRef1;
        undef $readRef2;
        undef $readCons1;
        undef $readCons2;

        $k += 1;    
      }
      $k = 1;
    }  
  }
}
close OUT1;
close OUT2;


