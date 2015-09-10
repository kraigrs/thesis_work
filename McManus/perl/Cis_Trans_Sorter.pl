#!/usr/bin/perl
# Joel McManus, 01/18/2010
# This script takes the tab-delimited data (see Table S2) and sorts the genes into divergence categories based on the q-values from the 
# binomial exact and Fisher exact test results.  This script sorts at 20, 50, 100, 200, and 500 read-hit thresholds.  Output includes a
# tab-delimited file for each sorting criteria (with the original data and the category for each gene) and a summary file of all results
# from that particular sorting.

use FileHandle;
use strict;

my $input;
my $qval_cutoff = 0.005;
my $readhit_cutoff = 20;

$input = "Cis_Trans_Stat_Table_Jan.txt";

my $OUTPUT = new FileHandle ">Cis_reg_"."$qval_cutoff"."_$readhit_cutoff".".txt" or die "can't open Cons_reg.txt";   #
$OUTPUT->printf("Gene\tHyb_M_Hits\tHyb_M_Hits_Adj\tHyb_S_Hits\tHyb_S_Hits_Adj\tHyb_Bin_pvalues\tHyb_Bin_qvalues\tPar_M_Hits\tPar_S_Hits\tPar_S_Hits_Adj\tPar_M_Hits_Norm\tPar_M_Hits_Adj\tPar_Bin_pvalues\tPar_Bin_qvalues\tLog2_M/S_Hyb\tLog2_M/S_Par\tFET_pvalues\tFET_qvalues\tHyb_Total\tCategory\n");
my $OUTPUT2 = new FileHandle ">Div_Cat_Count_"."$qval_cutoff".".txt" or die "can't open ouput 2";
$OUTPUT2->printf("$qval_cutoff=qvalue\t$readhit_cutoff=expression_level\nCIS\tTRANS\tCIS_PLUS_TRANS\tCIS_X_TRANS\tCOMPENSATORY\tCONSERVED\tCLASSIFIED\tPASS_COUNT\n");

my @f;
my $Mix_total;
my $line;
my $Hyb_M_hits;
my $Hyb_S_hits;
my $Hyb_qvalues;
my $Mix_M_hits;
my $Mix_S_hits;
my $Mix_qvalues;
my $FET_qvalues;
my $LOG2_P_Mel_Sec;
my $LOG2_H_Mel_Sec;
my $Mix_total;
my $cis;
my $gene;

my %CIS;
my %TRANS;
my %C_PLUS_T;
my %CIS_X_TRANS;
my %COMPENSATORY;
my %CONSERVED;
my $PASS_COUNT = 0;

open(RAWDATA, "$input") or die("couldn't open $input"); #tab-delimited post Excel file (mac line endings).  
## This will be the all_genes_input.txt file.  These are all the genes with 20+ total read-hits from parental sample.
<RAWDATA>; ### Should eat the header line ####
while ($line =<RAWDATA>) {
	chomp($line);
	@f = split(/\t/,$line);
	$gene = $f[0];
	$Mix_M_hits = $f[11];
	$Mix_S_hits = $f[9];
	$Mix_total = $Mix_M_hits + $Mix_S_hits;
	if ($Mix_total >= $readhit_cutoff) {
		$PASS_COUNT++;
	}
	$Mix_qvalues = $f[13];
	$Hyb_qvalues = $f[6];
	$FET_qvalues = $f[-2];
	$LOG2_P_Mel_Sec = $f[-4];
	$LOG2_H_Mel_Sec = $f[-5];
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS\n");
    	$CIS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tTRANS\n");
    	$TRANS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
		if ($LOG2_P_Mel_Sec < 0) {
			if ($LOG2_H_Mel_Sec < 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	} elsif ($LOG2_P_Mel_Sec > 0) {
    		if ($LOG2_H_Mel_Sec > 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	}
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS_X_TRANS_COMP\n");
    	$COMPENSATORY{$gene} = $line;
	}
		if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCONSERVED\n");
    	$CONSERVED{$gene} = $line;
	}
}
my $ciscount = keys %CIS;
my $transcount = keys %TRANS;
my $cisplustrans = keys %C_PLUS_T;
my $cisbytrans = keys %CIS_X_TRANS;
my $compcount = keys %COMPENSATORY;
my $conserved = keys %CONSERVED;
my $classified = $ciscount + $transcount + $cisplustrans + $cisbytrans + $compcount + $conserved;

my $cispct = sprintf("%.3f", (($ciscount/$PASS_COUNT)*100));
printf("$cispct\t");
my $transpct = sprintf("%.3f", (($transcount/$PASS_COUNT)*100));
printf("$transpct\t");
my $cisplustranspct = sprintf("%.3f", (($cisplustrans/$PASS_COUNT)*100));
printf("$cisplustranspct\t");
my $cisbytranspct = sprintf("%.3f", (($cisbytrans/$PASS_COUNT)*100));
printf("$cisbytranspct\t");
my $comppct = sprintf("%.3f", (($compcount/$PASS_COUNT)*100));
printf("$comppct\t");
my $conservedpct = sprintf("%.3f", (($conserved/$PASS_COUNT)*100));
printf("$conservedpct\t");
my $classifiedpct = sprintf("%.3f", (($classified/$PASS_COUNT)*100));
printf("$classifiedpct\n");

printf("$PASS_COUNT genes surpassed $readhit_cutoff reads\n");
$OUTPUT2->printf("$ciscount\t$transcount\t$cisplustrans\t$cisbytrans\t$compcount\t$conserved\t$classified\t$PASS_COUNT\n\n");
#$OUTPUT2->printf("$cispct\t$transpct\t$cisplustranspct\t$cisbytranspct\t$$comppct\t$conservedpct\t$classifiedpct\t$PASS_COUNT\n");
my $readhit_cutoff = 50;

$input = "Cis_Trans_Stat_Table_Jan_gt50.txt";

my $OUTPUT = new FileHandle ">Cis_reg_"."$qval_cutoff"."_$readhit_cutoff".".txt" or die "can't open Cons_reg.txt";   #
$OUTPUT->printf("Gene\tHyb_M_Hits\tHyb_M_Hits_Adj\tHyb_S_Hits\tHyb_S_Hits_Adj\tHyb_Bin_pvalues\tHyb_Bin_qvalues\tPar_M_Hits\tPar_S_Hits\tPar_S_Hits_Adj\tPar_M_Hits_Norm\tPar_M_Hits_Adj\tPar_Bin_pvalues\tPar_Bin_qvalues\tLog2_M/S_Hyb\tLog2_M/S_Par\tFET_pvalues\tFET_qvalues\tHyb_Total\tCategory\n");
$OUTPUT2->printf("$qval_cutoff=qvalue\t$readhit_cutoff=expression_level\nCIS\tTRANS\tCIS_PLUS_TRANS\tCIS_X_TRANS\tCOMPENSATORY\tCONSERVED\tCLASSIFIED\tPASS_COUNT\n");

my @f;
my $Mix_total;
my $line;
my $Hyb_M_hits;
my $Hyb_S_hits;
my $Hyb_qvalues;
my $Mix_M_hits;
my $Mix_S_hits;
my $Mix_qvalues;
my $FET_qvalues;
my $LOG2_P_Mel_Sec;
my $LOG2_H_Mel_Sec;
my $Mix_total;
my $cis;
my $gene;

my %CIS;
my %TRANS;
my %C_PLUS_T;
my %CIS_X_TRANS;
my %COMPENSATORY;
my %CONSERVED;
my $PASS_COUNT = 0;

open(RAWDATA, "$input") or die("couldn't open $input"); #tab-delimited post Excel file (mac line endings).  
## This will be the all_genes_input.txt file.  These are all the genes with 20+ total read-hits from parental sample.
<RAWDATA>; ### Should eat the header line ####
while ($line =<RAWDATA>) {
	chomp($line);
	@f = split(/\t/,$line);
	$gene = $f[0];
	$Mix_M_hits = $f[11];
	$Mix_S_hits = $f[9];
	$Mix_total = $Mix_M_hits + $Mix_S_hits;
	if ($Mix_total >= $readhit_cutoff) {
		$PASS_COUNT++;
	}
	$Mix_qvalues = $f[13];
	$Hyb_qvalues = $f[6];
	$FET_qvalues = $f[-2];
	$LOG2_P_Mel_Sec = $f[-4];
	$LOG2_H_Mel_Sec = $f[-5];
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS\n");
    	$CIS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tTRANS\n");
    	$TRANS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
		if ($LOG2_P_Mel_Sec < 0) {
			if ($LOG2_H_Mel_Sec < 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	} elsif ($LOG2_P_Mel_Sec > 0) {
    		if ($LOG2_H_Mel_Sec > 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	}
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS_X_TRANS_COMP\n");
    	$COMPENSATORY{$gene} = $line;
	}
		if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCONSERVED\n");
    	$CONSERVED{$gene} = $line;
	}
}
my $ciscount = keys %CIS;
my $transcount = keys %TRANS;
my $cisplustrans = keys %C_PLUS_T;
my $cisbytrans = keys %CIS_X_TRANS;
my $compcount = keys %COMPENSATORY;
my $conserved = keys %CONSERVED;
my $classified = $ciscount + $transcount + $cisplustrans + $cisbytrans + $compcount + $conserved;

my $cispct = sprintf("%.3f", (($ciscount/$PASS_COUNT)*100));
printf("$cispct\t");
my $transpct = sprintf("%.3f", (($transcount/$PASS_COUNT)*100));
printf("$transpct\t");
my $cisplustranspct = sprintf("%.3f", (($cisplustrans/$PASS_COUNT)*100));
printf("$cisplustranspct\t");
my $cisbytranspct = sprintf("%.3f", (($cisbytrans/$PASS_COUNT)*100));
printf("$cisbytranspct\t");
my $comppct = sprintf("%.3f", (($compcount/$PASS_COUNT)*100));
printf("$comppct\t");
my $conservedpct = sprintf("%.3f", (($conserved/$PASS_COUNT)*100));
printf("$conservedpct\t");
my $classifiedpct = sprintf("%.3f", (($classified/$PASS_COUNT)*100));
printf("$classifiedpct\n");

printf("$PASS_COUNT genes surpassed $readhit_cutoff reads\n");
$OUTPUT2->printf("$ciscount\t$transcount\t$cisplustrans\t$cisbytrans\t$compcount\t$conserved\t$classified\t$PASS_COUNT\n\n");
#$OUTPUT2->printf("$cispct\t$transpct\t$cisplustranspct\t$cisbytranspct\t$$comppct\t$conservedpct\t$classifiedpct\t$PASS_COUNT\n");
my $readhit_cutoff = 100;

$input = "Cis_Trans_Stat_Table_Jan_gt100.txt";

my $OUTPUT = new FileHandle ">Cis_reg_"."$qval_cutoff"."_$readhit_cutoff".".txt" or die "can't open Cons_reg.txt";   #
$OUTPUT->printf("Gene\tHyb_M_Hits\tHyb_M_Hits_Adj\tHyb_S_Hits\tHyb_S_Hits_Adj\tHyb_Bin_pvalues\tHyb_Bin_qvalues\tPar_M_Hits\tPar_S_Hits\tPar_S_Hits_Adj\tPar_M_Hits_Norm\tPar_M_Hits_Adj\tPar_Bin_pvalues\tPar_Bin_qvalues\tLog2_M/S_Hyb\tLog2_M/S_Par\tFET_pvalues\tFET_qvalues\tHyb_Total\tCategory\n");
$OUTPUT2->printf("$qval_cutoff=qvalue\t$readhit_cutoff=expression_level\nCIS\tTRANS\tCIS_PLUS_TRANS\tCIS_X_TRANS\tCOMPENSATORY\tCONSERVED\tCLASSIFIED\tPASS_COUNT\n");

my @f;
my $Mix_total;
my $line;
my $Hyb_M_hits;
my $Hyb_S_hits;
my $Hyb_qvalues;
my $Mix_M_hits;
my $Mix_S_hits;
my $Mix_qvalues;
my $FET_qvalues;
my $LOG2_P_Mel_Sec;
my $LOG2_H_Mel_Sec;
my $Mix_total;
my $cis;
my $gene;

my %CIS;
my %TRANS;
my %C_PLUS_T;
my %CIS_X_TRANS;
my %COMPENSATORY;
my %CONSERVED;
my $PASS_COUNT = 0;

open(RAWDATA, "$input") or die("couldn't open $input"); #tab-delimited post Excel file (mac line endings).  
## This will be the all_genes_input.txt file.  These are all the genes with 20+ total read-hits from parental sample.
<RAWDATA>; ### Should eat the header line ####
while ($line =<RAWDATA>) {
	chomp($line);
	@f = split(/\t/,$line);
	$gene = $f[0];
	$Mix_M_hits = $f[11];
	$Mix_S_hits = $f[9];
	$Mix_total = $Mix_M_hits + $Mix_S_hits;
	if ($Mix_total >= $readhit_cutoff) {
		$PASS_COUNT++;
	}
	$Mix_qvalues = $f[13];
	$Hyb_qvalues = $f[6];
	$FET_qvalues = $f[-2];
	$LOG2_P_Mel_Sec = $f[-4];
	$LOG2_H_Mel_Sec = $f[-5];
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS\n");
    	$CIS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tTRANS\n");
    	$TRANS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
		if ($LOG2_P_Mel_Sec < 0) {
			if ($LOG2_H_Mel_Sec < 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	} elsif ($LOG2_P_Mel_Sec > 0) {
    		if ($LOG2_H_Mel_Sec > 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	}
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS_X_TRANS_COMP\n");
    	$COMPENSATORY{$gene} = $line;
	}
		if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCONSERVED\n");
    	$CONSERVED{$gene} = $line;
	}
}
my $ciscount = keys %CIS;
my $transcount = keys %TRANS;
my $cisplustrans = keys %C_PLUS_T;
my $cisbytrans = keys %CIS_X_TRANS;
my $compcount = keys %COMPENSATORY;
my $conserved = keys %CONSERVED;
my $classified = $ciscount + $transcount + $cisplustrans + $cisbytrans + $compcount + $conserved;

my $cispct = sprintf("%.3f", (($ciscount/$PASS_COUNT)*100));
printf("$cispct\t");
my $transpct = sprintf("%.3f", (($transcount/$PASS_COUNT)*100));
printf("$transpct\t");
my $cisplustranspct = sprintf("%.3f", (($cisplustrans/$PASS_COUNT)*100));
printf("$cisplustranspct\t");
my $cisbytranspct = sprintf("%.3f", (($cisbytrans/$PASS_COUNT)*100));
printf("$cisbytranspct\t");
my $comppct = sprintf("%.3f", (($compcount/$PASS_COUNT)*100));
printf("$comppct\t");
my $conservedpct = sprintf("%.3f", (($conserved/$PASS_COUNT)*100));
printf("$conservedpct\t");
my $classifiedpct = sprintf("%.3f", (($classified/$PASS_COUNT)*100));
printf("$classifiedpct\n");

printf("$PASS_COUNT genes surpassed $readhit_cutoff reads\n");
$OUTPUT2->printf("$ciscount\t$transcount\t$cisplustrans\t$cisbytrans\t$compcount\t$conserved\t$classified\t$PASS_COUNT\n\n");
#$OUTPUT2->printf("$cispct\t$transpct\t$cisplustranspct\t$cisbytranspct\t$$comppct\t$conservedpct\t$classifiedpct\t$PASS_COUNT\n");

my $readhit_cutoff = 200;

$input = "Cis_Trans_Stat_Table_Jan_gt200.txt";

my $OUTPUT = new FileHandle ">Cis_reg_"."$qval_cutoff"."_$readhit_cutoff".".txt" or die "can't open Cons_reg.txt";   #
$OUTPUT->printf("Gene\tHyb_M_Hits\tHyb_M_Hits_Adj\tHyb_S_Hits\tHyb_S_Hits_Adj\tHyb_Bin_pvalues\tHyb_Bin_qvalues\tPar_M_Hits\tPar_S_Hits\tPar_S_Hits_Adj\tPar_M_Hits_Norm\tPar_M_Hits_Adj\tPar_Bin_pvalues\tPar_Bin_qvalues\tLog2_M/S_Hyb\tLog2_M/S_Par\tFET_pvalues\tFET_qvalues\tHyb_Total\tCategory\n");
$OUTPUT2->printf("$qval_cutoff=qvalue\t$readhit_cutoff=expression_level\nCIS\tTRANS\tCIS_PLUS_TRANS\tCIS_X_TRANS\tCOMPENSATORY\tCONSERVED\tCLASSIFIED\tPASS_COUNT\n");

my @f;
my $Mix_total;
my $line;
my $Hyb_M_hits;
my $Hyb_S_hits;
my $Hyb_qvalues;
my $Mix_M_hits;
my $Mix_S_hits;
my $Mix_qvalues;
my $FET_qvalues;
my $LOG2_P_Mel_Sec;
my $LOG2_H_Mel_Sec;
my $Mix_total;
my $cis;
my $gene;

my %CIS;
my %TRANS;
my %C_PLUS_T;
my %CIS_X_TRANS;
my %COMPENSATORY;
my %CONSERVED;
my $PASS_COUNT = 0;

open(RAWDATA, "$input") or die("couldn't open $input"); #tab-delimited post Excel file (mac line endings).  
## This will be the all_genes_input.txt file.  These are all the genes with 20+ total read-hits from parental sample.
<RAWDATA>; ### Should eat the header line ####
while ($line =<RAWDATA>) {
	chomp($line);
	@f = split(/\t/,$line);
	$gene = $f[0];
	$Mix_M_hits = $f[11];
	$Mix_S_hits = $f[9];
	$Mix_total = $Mix_M_hits + $Mix_S_hits;
	if ($Mix_total >= $readhit_cutoff) {
		$PASS_COUNT++;
	}
	$Mix_qvalues = $f[13];
	$Hyb_qvalues = $f[6];
	$FET_qvalues = $f[-2];
	$LOG2_P_Mel_Sec = $f[-4];
	$LOG2_H_Mel_Sec = $f[-5];
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS\n");
    	$CIS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tTRANS\n");
    	$TRANS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
		if ($LOG2_P_Mel_Sec < 0) {
			if ($LOG2_H_Mel_Sec < 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	} elsif ($LOG2_P_Mel_Sec > 0) {
    		if ($LOG2_H_Mel_Sec > 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	}
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS_X_TRANS_COMP\n");
    	$COMPENSATORY{$gene} = $line;
	}
		if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCONSERVED\n");
    	$CONSERVED{$gene} = $line;
	}
}
my $ciscount = keys %CIS;
my $transcount = keys %TRANS;
my $cisplustrans = keys %C_PLUS_T;
my $cisbytrans = keys %CIS_X_TRANS;
my $compcount = keys %COMPENSATORY;
my $conserved = keys %CONSERVED;
my $classified = $ciscount + $transcount + $cisplustrans + $cisbytrans + $compcount + $conserved;

my $cispct = sprintf("%.3f", (($ciscount/$PASS_COUNT)*100));
printf("$cispct\t");
my $transpct = sprintf("%.3f", (($transcount/$PASS_COUNT)*100));
printf("$transpct\t");
my $cisplustranspct = sprintf("%.3f", (($cisplustrans/$PASS_COUNT)*100));
printf("$cisplustranspct\t");
my $cisbytranspct = sprintf("%.3f", (($cisbytrans/$PASS_COUNT)*100));
printf("$cisbytranspct\t");
my $comppct = sprintf("%.3f", (($compcount/$PASS_COUNT)*100));
printf("$comppct\t");
my $conservedpct = sprintf("%.3f", (($conserved/$PASS_COUNT)*100));
printf("$conservedpct\t");
my $classifiedpct = sprintf("%.3f", (($classified/$PASS_COUNT)*100));
printf("$classifiedpct\n");

printf("$PASS_COUNT genes surpassed $readhit_cutoff reads\n");
$OUTPUT2->printf("$ciscount\t$transcount\t$cisplustrans\t$cisbytrans\t$compcount\t$conserved\t$classified\t$PASS_COUNT\n\n");
#$OUTPUT2->printf("$cispct\t$transpct\t$cisplustranspct\t$cisbytranspct\t$$comppct\t$conservedpct\t$classifiedpct\t$PASS_COUNT\n");
my $readhit_cutoff = 500;

$input = "Cis_Trans_Stat_Table_Jan_gt200.txt";

my $OUTPUT = new FileHandle ">Cis_reg_"."$qval_cutoff"."_$readhit_cutoff".".txt" or die "can't open Cons_reg.txt";   #
$OUTPUT->printf("Gene\tHyb_M_Hits\tHyb_M_Hits_Adj\tHyb_S_Hits\tHyb_S_Hits_Adj\tHyb_Bin_pvalues\tHyb_Bin_qvalues\tPar_M_Hits\tPar_S_Hits\tPar_S_Hits_Adj\tPar_M_Hits_Norm\tPar_M_Hits_Adj\tPar_Bin_pvalues\tPar_Bin_qvalues\tLog2_M/S_Hyb\tLog2_M/S_Par\tFET_pvalues\tFET_qvalues\tHyb_Total\tCategory\n");
$OUTPUT2->printf("$qval_cutoff=qvalue\t$readhit_cutoff=expression_level\nCIS\tTRANS\tCIS_PLUS_TRANS\tCIS_X_TRANS\tCOMPENSATORY\tCONSERVED\tCLASSIFIED\tPASS_COUNT\n");

my @f;
my $Mix_total;
my $line;
my $Hyb_M_hits;
my $Hyb_S_hits;
my $Hyb_qvalues;
my $Mix_M_hits;
my $Mix_S_hits;
my $Mix_qvalues;
my $FET_qvalues;
my $LOG2_P_Mel_Sec;
my $LOG2_H_Mel_Sec;
my $Mix_total;
my $cis;
my $gene;

my %CIS;
my %TRANS;
my %C_PLUS_T;
my %CIS_X_TRANS;
my %COMPENSATORY;
my %CONSERVED;
my $PASS_COUNT = 0;

open(RAWDATA, "$input") or die("couldn't open $input"); #tab-delimited post Excel file (mac line endings).  
## This will be the all_genes_input.txt file.  These are all the genes with 20+ total read-hits from parental sample.
<RAWDATA>; ### Should eat the header line ####
while ($line =<RAWDATA>) {
	chomp($line);
	@f = split(/\t/,$line);
	$gene = $f[0];
	$Mix_M_hits = $f[11];
	$Mix_S_hits = $f[9];
	$Mix_total = $Mix_M_hits + $Mix_S_hits;
	if ($Mix_total >= $readhit_cutoff) {
		$PASS_COUNT++;
	}
	$Mix_qvalues = $f[13];
	$Hyb_qvalues = $f[6];
	$FET_qvalues = $f[-2];
	$LOG2_P_Mel_Sec = $f[-4];
	$LOG2_H_Mel_Sec = $f[-5];
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS\n");
    	$CIS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tTRANS\n");
    	$TRANS{$gene} = $line;
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues < $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
		if ($LOG2_P_Mel_Sec < 0) {
			if ($LOG2_H_Mel_Sec < 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	} elsif ($LOG2_P_Mel_Sec > 0) {
    		if ($LOG2_H_Mel_Sec > 0) {
    			$OUTPUT->printf("$line\tC_PLUS_T\n");
    			$C_PLUS_T{$gene} = $line;
    		} else {
    			$OUTPUT->printf("$line\tC_X_T\n");
    			$CIS_X_TRANS{$gene} = $line;
    		}
    	}
	}
	if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues < $qval_cutoff) && ($Hyb_qvalues < $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCIS_X_TRANS_COMP\n");
    	$COMPENSATORY{$gene} = $line;
	}
		if (($Mix_total >= $readhit_cutoff) && ($Mix_qvalues > $qval_cutoff) && ($FET_qvalues > $qval_cutoff) && ($Hyb_qvalues > $qval_cutoff)) {
    	$OUTPUT->printf("$line\tCONSERVED\n");
    	$CONSERVED{$gene} = $line;
	}
}
my $ciscount = keys %CIS;
my $transcount = keys %TRANS;
my $cisplustrans = keys %C_PLUS_T;
my $cisbytrans = keys %CIS_X_TRANS;
my $compcount = keys %COMPENSATORY;
my $conserved = keys %CONSERVED;
my $classified = $ciscount + $transcount + $cisplustrans + $cisbytrans + $compcount + $conserved;

my $cispct = sprintf("%.3f", (($ciscount/$PASS_COUNT)*100));
printf("$cispct\t");
my $transpct = sprintf("%.3f", (($transcount/$PASS_COUNT)*100));
printf("$transpct\t");
my $cisplustranspct = sprintf("%.3f", (($cisplustrans/$PASS_COUNT)*100));
printf("$cisplustranspct\t");
my $cisbytranspct = sprintf("%.3f", (($cisbytrans/$PASS_COUNT)*100));
printf("$cisbytranspct\t");
my $comppct = sprintf("%.3f", (($compcount/$PASS_COUNT)*100));
printf("$comppct\t");
my $conservedpct = sprintf("%.3f", (($conserved/$PASS_COUNT)*100));
printf("$conservedpct\t");
my $classifiedpct = sprintf("%.3f", (($classified/$PASS_COUNT)*100));
printf("$classifiedpct\n");

printf("$PASS_COUNT genes surpassed $readhit_cutoff reads\n");
$OUTPUT2->printf("$ciscount\t$transcount\t$cisplustrans\t$cisbytrans\t$compcount\t$conserved\t$classified\t$PASS_COUNT\n\n");