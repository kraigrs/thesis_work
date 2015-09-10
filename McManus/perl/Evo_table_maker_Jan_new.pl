#!/usr/bin/perl
#

#---------SCRIPT PARAMETERS
use POSIX qw(ceil floor);
use strict;
use FileHandle;

#---------VARIABLE DEFINITIONS

my $READS1 = "Hyb_mel_nogap_unq_cons_gene_hit_counts_Jan_new.txt";
my $READS2 = "Hyb_sec_lift_nogap_unq_cons_gene_hit_counts_Jan_new.txt";
my $READS3 = "Mix_mel_nogap_unq_cons_gene_hit_counts_Jan_new.txt";
my $READS4 = "Mix_sec_lift_nogap_unq_cons_gene_hit_counts_Jan_new.txt";

my $CIS_TRANS_TABLE1 = new FileHandle ">Cis_Trans_Evo_Table_Jan_new.txt" or die "can't open OUTFILE";
$CIS_TRANS_TABLE1->printf("Gene\tHyb_M_Hits\tHyb_S_Hits\tPar_M_Hits\tPar_S_Hits\tPar_M_Hits_Norm\n");
my $CIS_TRANS_TABLE2 = new FileHandle ">Cis_Trans_Evo_ChrM_Table_Jan_new.txt" or die "can't open OUTFILE2";
$CIS_TRANS_TABLE2->printf("Gene\tHyb_M_Hits\tHyb_S_Hits\tPar_M_Hits\tPar_S_Hits\n");


my %HYB_MEL_COUNT;
my %HYB_SEC_COUNT;
my %MIX_MEL_COUNT;
my %MIX_SEC_COUNT;
my %HYB_MEL_CHRM_COUNT;
my %HYB_SEC_CHRM_COUNT;
my %MIX_MEL_CHRM_COUNT;
my %MIX_SEC_CHRM_COUNT;

my %ALL_GENES;
my %CHRM_GENES = ("CG34063", 1, "CG34067", 1, "CG34069", 1, "CG34072", 1, "CG34073", 1, "CG34074", 1, "CG34076", 1, "CG34083", 1, "CG34085", 1, "CG34086", 1, "CG34089", 1, "CG34090", 1, "CG34092", 1);


my $gene;
my $Hyb_Mel_Total = 0;
my $Hyb_Sec_Total = 0;
my $Mix_Mel_Total = 0;
my $Mix_Sec_Total = 0;
my $Hyb_Mel_ChrM_Total = 0;
my $Hyb_Sec_ChrM_Total = 0;
my $Mix_Mel_ChrM_Total = 0;
my $Mix_Sec_ChrM_Total = 0;

my $MIX_MEL_NORMAL;
my $MIX_MEL_NORMAL_UP;
my $LOG2_H_MEL_SEC;
my $LOG2_P_MEL_SEC;


open(READS1, "<$READS1") or die("can't open $READS1");

my @f;
my $line;

while($line= <READS1>){
  chomp($line);
  @f = split /\s+/,$line;
  $gene = $f[0];  
  if (exists $CHRM_GENES{$gene}) {### filters out mitochondrial genes into a new file
  	$HYB_MEL_CHRM_COUNT{$gene} = $f[1];
  	$Hyb_Mel_ChrM_Total +=$f[1];
  } else {
  	$HYB_MEL_COUNT{$gene} = $f[1];
  	$ALL_GENES{$gene} = 1;
  	$Hyb_Mel_Total +=$f[1];
  }
}
close (READS1);
print "$Hyb_Mel_Total Hybrid Melanogaster Gene hits\n";
open(READS2, "<$READS2") or die("can't open $READS2");

my @f;
my $line;

while($line= <READS2>){
  chomp($line);
  @f = split /\s+/,$line;
  $gene = $f[0];  
  if (exists $CHRM_GENES{$gene}) {### filters out mitochondrial genes into a new file
  	$HYB_SEC_CHRM_COUNT{$gene} = $f[1];
  	$Hyb_Sec_ChrM_Total +=$f[1];
  } else {
  	$HYB_SEC_COUNT{$gene} = $f[1];
  	$ALL_GENES{$gene} = 1;
  	$Hyb_Sec_Total +=$f[1];
  }
}
close (READS2);
print "$Hyb_Sec_Total Hybrid Sechellia Gene hits\n";
open(READS3, "<$READS3") or die("can't open $READS3");
my @f;
my $line;

while($line= <READS3>){
  chomp($line);
  @f = split /\s+/,$line;
  $gene = $f[0];  
  if (exists $CHRM_GENES{$gene}) {### filters out mitochondrial genes into a new file
  	$MIX_MEL_CHRM_COUNT{$gene} = $f[1];
  	$Mix_Mel_ChrM_Total +=$f[1];
  } else {
  	$MIX_MEL_COUNT{$gene} = $f[1];
  	$ALL_GENES{$gene} = 1;
  	$Mix_Mel_Total +=$f[1];
  }
}
close (READS3);
print "$Mix_Mel_Total Parental Melanogaster Gene hits\n";
open(READS4, "<$READS4") or die("can't open $READS4");
my @f;
my $line;

while($line= <READS4>){
  chomp($line);
  @f = split /\s+/,$line;
  $gene = $f[0];  
  if (exists $CHRM_GENES{$gene}) {### filters out mitochondrial genes into a new file
  	$MIX_SEC_CHRM_COUNT{$gene} = $f[1];
  	$Mix_Sec_ChrM_Total +=$f[1];
  } else {
  	$MIX_SEC_COUNT{$gene} = $f[1];
  	$ALL_GENES{$gene} = 1;
  	$Mix_Sec_Total +=$f[1];
  }
}
close (READS4);
print "$Mix_Sec_Total Parental Sechellia Gene hits\n";

my $Hyb_Mel_Pct_chrM = ($Hyb_Mel_ChrM_Total / ($Hyb_Mel_ChrM_Total + $Hyb_Sec_ChrM_Total)) * 100;
my $Mix_Mel_Pct_chrM = ($Mix_Mel_ChrM_Total / ($Mix_Mel_ChrM_Total + $Mix_Sec_ChrM_Total)) * 100;
my $Hyb_Mel_over_Sec = $Hyb_Mel_Total / $Hyb_Sec_Total;
my $Mix_Mel_over_Sec = $Mix_Mel_Total / $Mix_Sec_Total;
my $Normal_factor = $Mix_Mel_over_Sec / $Hyb_Mel_over_Sec;

print "$Hyb_Mel_Pct_chrM percent of Hybrid mitochondrial reads from Melanogaster\n";
print "$Mix_Mel_Pct_chrM percent of Mixed Parental mitochondrial reads from Melanogaster\n";

foreach my $gene (keys %ALL_GENES) {
	$MIX_MEL_NORMAL = $MIX_MEL_COUNT{$gene} / $Normal_factor;
	$MIX_MEL_NORMAL_UP = ceil($MIX_MEL_NORMAL);
	#$LOG2_H_MEL_SEC = (log(($HYB_MEL_COUNT{$gene}/$HYB_SEC_COUNT{$gene}))/log(2));
	#$LOG2_P_MEL_SEC = (log(($MIX_MEL_NORMAL_UP/$MIX_SEC_COUNT{$gene}))/log(2));
	if (($MIX_MEL_NORMAL_UP + $MIX_SEC_COUNT{$gene}) >= 20) {
	$CIS_TRANS_TABLE1->printf("$gene\t$HYB_MEL_COUNT{$gene}\t$HYB_SEC_COUNT{$gene}\t$MIX_MEL_COUNT{$gene}\t$MIX_SEC_COUNT{$gene}\t$MIX_MEL_NORMAL_UP\n");
	}
}

foreach my $gene (keys %CHRM_GENES) {
	$CIS_TRANS_TABLE2->printf("$gene\t$HYB_MEL_CHRM_COUNT{$gene}\t$HYB_SEC_CHRM_COUNT{$gene}\t$MIX_MEL_CHRM_COUNT{$gene}\t$MIX_SEC_CHRM_COUNT{$gene}\n");
}



