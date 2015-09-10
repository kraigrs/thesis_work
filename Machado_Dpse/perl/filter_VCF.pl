#!/usr/bin/perl

######################################################################################
# 
# 11/06/2013
#
# filter_VCF.pl
#
# Purpose: read in pileup from SAMtools to call SNPs (or mutations)
# 
# Input: SAMtools alignment pileup and VCF file         
#
# Output: based on certain criteria, return highest confidence variants
#
######################################################################################

use strict;
use warnings;

my$pileup = $ARGV[0];
my$VCF = $ARGV[1];
my$output = $ARGV[2];

my$line; my$locus; my$pos; my$covseqs; my$call; my$match; my$As; my$Cs; my$Gs; my$Ts; 
my$phred; my$sum; my$ref; my$alt; my$insertions; my$deletions; my$Q; my$mismatch; my$avg;
my$VCF_snp_counter; my$pile_snp_counter; my$VCF_indel_counter; my$pile_indel_counter; 
my$common_snps; my$common_indels; my$i; my$temp;
my@elements; my@bases; my@quals;
my%confident;

sub prob
{
  my$num = shift;
  my$exponent = -$num/10;
  my$result = 10**$exponent;
  return $result;
}

sub lower
{
  my$num = shift;
  my$result = 0.5*$num-sqrt(0.25*$num);
  return $result;
}

sub upper
{
  my$num = shift;
  my$result = 0.5*$num+sqrt(0.25*$num);
  return $result;
}

$pile_snp_counter = 0; 
$pile_indel_counter = 0; 
$VCF_snp_counter = 0; 
$VCF_indel_counter = 0; 

open(PILE,"$pileup") or die "\nError opening $pileup\n";
while(<PILE>)
{
  #$i += 1;
  #if($i % 1000000 == 0)
  #{
  #  print "Line: $i\tSNPs: $pile_snp_counter\tindels: $pile_indel_counter\n";
  #}

  chomp;
  $line = $_;
  unless(/^#/)
  {
    @elements = split("\t",$_);
    $locus = $elements[0];
    $pos = $elements[1];
    $covseqs = $elements[3];
    $call = $elements[4];
    $Q = $elements[5];
  
    $match = 0;
    $mismatch = 0;
    $insertions = 0;
    $deletions = 0;
    $sum = 0;
    $avg = 0;
    $As = 0;
    $Cs = 0;
    $Gs = 0;
    $Ts = 0; 
    
    @bases = split("",$call);
    foreach(@bases)
    {
      if(/[\.,]/){$match+=1;}
      elsif(/[\+]/){$insertions+=1;}
      elsif(/[\-]/){$deletions+=1;}
      elsif(/[Aa]/){$As+=1;}
      elsif(/[Cc]/){$Cs+=1;}
      elsif(/[Gg]/){$Gs+=1;}
      elsif(/[Tt]/){$Ts+=1;}
    }

    @quals = split("",$Q);
    foreach(@quals)
    {
      $phred = ord($_)-33;
      $sum += prob($phred);
    }

    $avg = $sum/$covseqs;
    #$mismatch = $As + $Cs + $Gs + $Ts;

    if($covseqs >= 20 && $avg < 0.01)
    {
      if($insertions + $deletions == 0)
      {
        if( ($As > 0 && $Cs+$Gs+$Ts == 0) ||
            ($Cs > 0 && $As+$Gs+$Ts == 0) ||
            ($Gs > 0 && $Cs+$As+$Ts == 0) ||
            ($Ts > 0 && $Cs+$Gs+$As == 0) )
        {
          #if( ($As>=lower($covseqs) && $As<=upper($covseqs)) || 
	  #    ($Cs>=lower($covseqs) && $Cs<=upper($covseqs)) || 
	  #    ($Gs>=lower($covseqs) && $Gs<=upper($covseqs)) || 
	  #    ($Ts>=lower($covseqs) && $Ts<=upper($covseqs)) )
	  #{
          #  $confident{$locus}{$pos} = "heterozygous_SNP";
          #  $pile_snp_counter += 1;
          #  #print "$line\n";
	  #}
          #elsif($As>upper($covseqs) || $Cs>upper($covseqs) || $Gs>upper($covseqs) || $Ts>upper($covseqs))
	  #{
          #  $confident{$locus}{$pos} = "homozygous_SNP";
          #  $pile_snp_counter += 1;
          #  #print "$line\n";
	  #}
          if($As==$covseqs || $Cs==$covseqs || $Gs==$covseqs || $Ts==$covseqs)
	  {
            $confident{$locus}{$pos} = "homozygous_SNP";
            $pile_snp_counter += 1;
            #print "$line\n";
	  }
        }
      }
      #elsif($insertions > 0 && $deletions == 0)
      #{
        #if($insertions >= lower($covseqs) && $insertions <= upper($covseqs))
	#{
        #  $confident{$locus}{$pos} = "heterozygous_insertion";
        #  $pile_indel_counter += 1;
        #  #print "$line\n";
	#}
        #elsif($insertions > upper($covseqs))
	#{
        #  $confident{$locus}{$pos} = "homozygous_insertion";
        #  $pile_indel_counter += 1;
        #  #print "$line\n";
	#}
        #if($insertions == $covseqs)
	#{
        #  $confident{$locus}{$pos} = "homozygous_insertion";
        #  $pile_indel_counter += 1;
        #  #print "$line\n";
	#}
      #}
      #elsif($insertions == 0 && $deletions > 0)
      #{
        #if($deletions >= lower($covseqs) && $deletions <= upper($covseqs))
	#{
        #  $confident{$locus}{$pos} = "heterozygous_deletion";
        #  $pile_indel_counter += 1;
        #  #print "$line\n";
	#}
        #elsif($deletions > upper($covseqs))
	#{
        #  $confident{$locus}{$pos} = "homozygous_deletion";
        #  $pile_indel_counter += 1;
        #  #print "$line\n";
	#}
        #if($deletions == $covseqs)
	#{
        #  $confident{$locus}{$pos} = "homozygous_deletion";
        #  $pile_indel_counter += 1;
          #print "$line\n";
	#}
      #}
    }
  }
}
close PILE;

open(OUT,"> $output") or die "Error writing to $output\n";

open(VCF,"$VCF") or die "\nError opening $VCF\n";
while(<VCF>)
{
  chomp;
  $line = $_;
  if(/^#/){print OUT "$line\n";}
  else
  {
    @elements = split("\t",$_);
    $locus = $elements[0];
    $pos = $elements[1];
    $ref = $elements[3];
    $alt = $elements[4];

    if(length($ref) > 1 || length($alt) > 1){$call = "indel"; $VCF_indel_counter += 1;}
    else{$call = "SNP"; $VCF_snp_counter += 1;}

    if($confident{$locus}{$pos})
    {
      print OUT "$line\n";
      if($confident{$locus}{$pos} =~ /SNP/){$common_snps += 1;}
      else{$common_indels += 1;}
    }
  }
}
close VCF;
close OUT;

print "\nTotals\n";

print "pileup SNPs: $pile_snp_counter\n";
print "common SNPs: $common_snps\n";
print "VCF SNPs: $VCF_snp_counter\n";

print "\npileup indels: $pile_indel_counter\n";
print "common indels: $common_indels\n";
print "VCF indels: $VCF_indel_counter\n";
