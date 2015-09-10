#!/usr/bin/perl
##  Joel McManus 07/23/2010
##  This script aligns RNA-seq data to find snp genomic and junction hits from Zhr and Z30.
###  !!!! NOTES !!!! Check the readIDs to make sure same lane, different date doesn't have same ID naming scheme!!!!!
use strict;
use FileHandle;
system "date";
printf("\n**********\n\tAligning reads to genome and splice junctions\n**********\n\n");
# This section aligns all the RNA-seq reads to the melanogaster zhr and z30 genomes and junctions. Reads that don't map to genome are written into "nomap" files and then aligned to junctions. Reads that don't map to genome or junction are then aligned to trans-spliced junctions.
my @readfiles = ("3144CAAXX_l6_r1_zhr_seq.txt", "3144CAAXX_l6_r2_zhr_seq.txt", "61FOVAAXX_l4_r2_z30_seq.txt", "61FOVAAXX_l4_r1_z30_seq.txt", "61F0CAAXX_l3_r2_z30xzhr_seq.txt", "61F0CAAXX_l3_r1_z30xzhr_seq.txt", "61F0CAAXX_l2_r2_zhrxz30_seq.txt", "61F0CAAXX_l2_r1_zhrxz30_seq.txt");
foreach (@readfiles) {
	my $nomap1 = $_; $nomap1 =~ s/seq.txt/nomap1.fq/;
	my $nomap2 = $_; $nomap2 =~ s/seq.txt/nomap2.fq/;
	my $nomap3 = $_; $nomap3 =~ s/seq.txt/nomap3.fq/;
	my $nomap4 = $_; $nomap4 =~ s/seq.txt/nomap4.fq/;
	my $map_zhr = $_; $map_zhr =~ s/seq.txt/zhr.dat/;
	my $map_z30 = $_; $map_z30 =~ s/seq.txt/z30.dat/;
	my $map_zhr2 = $_; $map_zhr2 =~ s/seq.txt/zhr2.dat/;
	my $map_z302 = $_; $map_z302 =~ s/seq.txt/z302.dat/;
	my $map_zhr3 = $_; $map_zhr3 =~ s/seq.txt/zhr3.dat/;
	my $map_z303 = $_; $map_z303 =~ s/seq.txt/z303.dat/;
	my $snp_zhr = $_; $snp_zhr =~ s/seq.txt/zhr_snp.dat/;
	my $snp_z30 = $_; $snp_z30 =~ s/seq.txt/z30_snp.dat/;
	my $snp_zhr2 = $_; $snp_zhr2 =~ s/seq.txt/zhr_snp2.dat/;
	my $snp_z302 = $_; $snp_z302 =~ s/seq.txt/z30_snp2.dat/;
	my $snp_zhr3 = $_; $snp_zhr3 =~ s/seq.txt/zhr_snp3.dat/;
	my $snp_z303 = $_; $snp_z303 =~ s/seq.txt/z30_snp3.dat/;
	my $mapjunc_zhr = $_; $mapjunc_zhr =~ s/seq.txt/junc_zhr.dat/;
	my $mapjunc_z30 = $_; $mapjunc_z30 =~ s/seq.txt/junc_z30.dat/;
	my $mapjunc_zhr2 = $_; $mapjunc_zhr2 =~ s/seq.txt/junc_zhr2.dat/;
	my $mapjunc_z302 = $_; $mapjunc_z302 =~ s/seq.txt/junc_z302.dat/;
	my $mapjunc_zhr3 = $_; $mapjunc_zhr3 =~ s/seq.txt/junc_zhr3.dat/;
	my $mapjunc_z303 = $_; $mapjunc_z303 =~ s/seq.txt/junc_z303.dat/;
	my $mapjunc_zhrz30_trans = $_; $mapjunc_zhrz30_trans =~ s/seq.txt/trans_junc.dat/;
	my $mapjunc_zhrz30_trans2 = $_; $mapjunc_zhrz30_trans2 =~ s/seq.txt/trans_junc2.dat/;
	my $mapjunc_zhrz30_trans3 = $_; $mapjunc_zhrz30_trans3 =~ s/seq.txt/trans_junc3.dat/;
	my $snpjunc_zhr = $_; $snpjunc_zhr =~ s/seq.txt/junc_zhr_snp.dat/;
	my $snpjunc_z30 = $_; $snpjunc_z30 =~ s/seq.txt/junc_z30_snp.dat/;
	my $snpjunc_zhr2 = $_; $snpjunc_zhr2 =~ s/seq.txt/junc_zhr_snp2.dat/;
	my $snpjunc_z302 = $_; $snpjunc_z302 =~ s/seq.txt/junc_z30_snp2.dat/;
	my $snpjunc_zhr3 = $_; $snpjunc_zhr3 =~ s/seq.txt/junc_zhr_snp3.dat/;
	my $snpjunc_z303 = $_; $snpjunc_z303 =~ s/seq.txt/junc_z30_snp3.dat/;
	my $trim = 76;
	my $mosaikreads;
	$mosaikreads = $_;
	printf("$_ and $mosaikreads\n\n");
	$mosaikreads =~ s/seq.txt/input.dat/;
	my $startfile;
	# Align and find snp reads to genomic and junction databases.
	while ($trim >= 37) {
		if ($trim == 76) {
			$startfile = $_;
		} else {
			$startfile = "$_"."_$trim";
		}
		printf("\n*********\nAligning reads in file $startfile to D. mel zhr and D. mel z30 genome databases.\n*********\n");
		system "MosaikBuild -q $startfile -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_zhr"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_AMB_13 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_zhr2"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_mix1.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_mix1_13 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_zhr3"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_mix2.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_mix2_13 -mhp 100 -act 20 -rur $nomap3";

		system "MosaikBuild -q $nomap3 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_z30"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_AMB_13 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_z302"."_$trim"." -ia ~/mosaik-aligner/targets/z30_mix1.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_mix1_13 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_z303"."_$trim"." -ia ~/mosaik-aligner/targets/z30_mix2.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_mix2_13 -mhp 100 -act 20 -rur $nomap3";

		system "MosaikBuild -q $startfile -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_z30"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_AMB_13 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_z302"."_$trim"." -ia ~/mosaik-aligner/targets/z30_mix1.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_mix1_13 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_z303"."_$trim"." -ia ~/mosaik-aligner/targets/z30_mix2.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_mix2_13 -mhp 100 -act 20 -rur $nomap3";

		system "MosaikBuild -q $nomap3 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_zhr"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_AMB_13 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_zhr2"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_mix1.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_mix1_13 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_zhr3"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_mix2.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_mix2_13 -mhp 100 -act 20 -rur $nomap3";

		### Align to junctions
		printf("\n*********\nAligning reads in file $startfile to D. mel zhr and D. mel z30 junction databases.\n*********\n");
		system "MosaikBuild -q $nomap3 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_zhr"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_zhr2"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB_mix_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_zhr3"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB_mix2_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap4";
		
		system "MosaikBuild -q $nomap4 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snpjunc_z30"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snpjunc_z302"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB_mix_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snpjunc_z303"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB_mix2_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap4";

		system "MosaikBuild -q $nomap3 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_z30"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_z302"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB_mix_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_z303"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB_mix2_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap4";

		system "MosaikBuild -q $nomap4 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snpjunc_zhr"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snpjunc_zhr2"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB_mix_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snpjunc_zhr3"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB_mix2_junc_l"."$trim".".dat -hs 13 -mm 0 -p 24 -mhp 100 -act 20 -rur $nomap4";
		### Align to trans-splicing junctions
		printf("\n*********\nAligning to trans-splice junctions for trim = $trim\n*********\n");
		system "MosaikBuild -q $nomap4 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_zhrz30_trans"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_z30_tj_"."$trim".".dat -hs 13 -mm 0 -p 24 -rur $nomap1";
		system "MosaikBuild -q $nomap1 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_zhrz30_trans2"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_z30_mix1_tj_"."$trim".".dat -hs 13 -mm 0 -p 24 -rur $nomap2";
		system "MosaikBuild -q $nomap2 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $mapjunc_zhrz30_trans3"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_z30_mix2_tj_"."$trim".".dat -hs 13 -mm 0 -p 24 -rur $nomap4";
		
		printf("\nFinished aligning to junctions for trim = $trim\n\n");
		$trim -= 13;
		if ($trim > 24) {
			printf("\n*********\nTrimming non-hitter reads to length=$trim for remapping.\n\n");
			system "perl fastq_trimmer.pl $nomap2 $trim > "."$_"."_$trim";	## trim reads for next round.
			}
		} ### Repeat above with trimmed reads and junctions.
} # align next reads file
system "date";