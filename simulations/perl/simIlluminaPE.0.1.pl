#!/usr/bin/perl

########################################################################
# 
# 05/11/2010
#
# simIlluminaPE0.1.pl
#
# Purpose: simulate sets of paired-end reads from reference databases
# 
# Input: a FASTA reference and SNP locations (with bp info)
#
# Output: millions of short sequence reads with various properties (mutation rates,
#         lengths, etc.)
# 
# Fixed: 05/11/2010 --> replaced last 3 loops with a single looping structure, more efficient  
#
#        05/12/2010 --> added both strand orientation simulated reads (basically doubling coverage)
#
# Notes: 05/11/2010 --> should consider replacing "find gene" loop with making a gene hash, 
#                       but may require too much memory
#      
#        05/12/2010 --> may be beneficial to also simulate reads other than those around SNPs.
#                       This may better represent reads that we later classify as "both".
#                       The simulation strategy would be something like this:
#                       1). Choose a gene (or chromosome) to start with, and generate reads starting 
#                           at the beginning of the reference sequence
#                       2). When a SNP is reached, start original process.
#
#        05/12/2010 --> also need to simulate reads in the opposite orientation!!!
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

my$ref = $ARGV[4];
my$out1 = $ARGV[5];
my$out2 = $ARGV[6];

my$seq; my$found = 0; my$i; my$j; my$k = 1;
my$reflength; my$posNeg; my$lengthTot; my$false1; my$false2;
my$readRef1; my$readCons1; my$readRef2; my$readCons2;
my$baseRef1; my$baseRef2; my$baseCons1; my$baseCons2;
my@nuc; my@nucNeg;

open(OUT1,">$out1") or die "\nError opening $out1\n";
open(OUT2,">$out2") or die "\nError opening $out2\n";

foreach $gene (sort keys %SNPhash)
{
  open(REF,"$ref") or die "\nError opening $ref\n";
  HERE: while(<REF>)    
  {
    chomp;
    if(/^>$gene\s+/)
    {      
      $found = 1; 
    }
    elsif(/^>.+/ && $found)
    {   
      last HERE;
    }
    elsif($found)
    {
      $seq = $seq . $_;  
    }
  }
  close REF;
  $found = 0;
  #print "$gene\n$seq\n";
  @nuc = split("",$seq);
  @nucNeg = reverse @nuc;
  undef $seq;
  $reflength = scalar@nuc;
  $lengthTot = 2*$length+$insert; # total length of paired-end read
  #print "Length of reference: $reflength\n";

  foreach $pos (sort keys %{$SNPhash{$gene}})
  {
    unless($lengthTot > $reflength)
    {
      # positive orientation

      for($i=max($pos-$lengthTot,0);$i<min($reflength-$lengthTot+1,$pos);$i++) # account for read exceeding reference
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
              
              $readRef1 = $readRef1 . $nuc[$i+$j+$length+$insert];
              if($i+$j+$length+$insert == $pos-1){$readCons1 = $readCons1 . $SNPhash{$gene}{$pos};}
              else{$readCons1 = $readCons1 . $nuc[$i+$j+$length+$insert];}
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
        print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\+\)\/1\n$readRef1\n";
        print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\+\)\/1\n$readCons1\n";

        print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\+\)\/2\n$readRef2\n";
        print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\+\)\/2\n$readCons2\n";

        undef $readRef1;
        undef $readRef2;
        undef $readCons1;
        undef $readCons2;

        $k += 1;    
      }
      $k = 1;
     
      # negative orientation

      $posNeg = $reflength - $pos + 1; # reverse the position of the SNP
      for($i=max($posNeg-$lengthTot,0);$i<min($reflength-$lengthTot+1,$posNeg);$i++) # account for read exceeding reference
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
              
              $readRef1 = $readRef1 . $nucNeg[$i+$j];
              if($i+$j == $posNeg-1){$readCons1 = $readCons1 . $SNPhash{$gene}{$pos};}
              else{$readCons1 = $readCons1 . $nucNeg[$i+$j];}
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
        print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\-\)\/1\n$readRef1\n";
        print OUT1 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\-\)\/1\n$readCons1\n";

        print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\-\)\/2\n$readRef2\n";
        print OUT2 "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\-\)\/2\n$readCons2\n";

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

# alternative gene hash structure

#open(REF,"$ref") or die "\nError opening $ref\n";
#while(<REF>)    
#{
#  chomp;
#  if(/^>(\S+)/)
#  {
#    $gene = $1;
#    undef $seq;
#  }
#  else
#  {
#    $seq = $seq . $_;
#    $geneHash{$gene} = $seq;  
#  }
#}
#close REF;

