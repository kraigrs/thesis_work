#!/bin/sh

##### lift exons for each comparison #####

# zhr and z30

#liftOver constitutive_regions_overlap_filtered2.bed dm3_zhr_AMBIG_all.chain constitutive_regions_overlap_filtered2.zhr.bed constitutive_regions_overlap_filtered2.zhr.un

#liftOver constitutive_regions_overlap_filtered2.bed dm3_z30_AMBIG_all.chain constitutive_regions_overlap_filtered2.z30.bed constitutive_regions_overlap_filtered2.z30.un

# sim

#liftOver constitutive_regions_overlap_filtered2.bed dm3ToDroSim1.over.chain constitutive_regions_overlap_filtered2.sim.bed constitutive_regions_overlap_filtered2.sim.un

#liftOver constitutive_regions_overlap_filtered2.sim.bed droSim1ToDroSim1_snp_indel.over.chain constitutive_regions_overlap_filtered2.sim_snp_indel.bed constitutive_regions_overlap_filtered2.sim_snp_indel.un

# sec

#liftOver constitutive_regions_overlap_filtered2.bed dm3ToDroSec1.over.chain constitutive_regions_overlap_filtered2.sec.bed constitutive_regions_overlap_filtered2.sec.un

#liftOver constitutive_regions_overlap_filtered2.sec.bed droSec1ToDroSec1_snp_indel.over.chain constitutive_regions_overlap_filtered2.sec_snp_indel.bed constitutive_regions_overlap_filtered2.sec_snp_indel.un

##### get sequences #####

# dm3

#perl get_sequences.pl constitutive_regions_overlap_filtered2.bed ../Dmel/dm3.fa constitutive_regions_overlap_filtered2.dm3.fa

# zhr and z30

#perl get_sequences.pl constitutive_regions_overlap_filtered2.zhr.bed zhr_AMBIG_all.fa constitutive_regions_overlap_filtered2.zhr.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2.z30.bed z30_AMBIG_all.fa constitutive_regions_overlap_filtered2.z30.fa

# sim

#perl get_sequences.pl constitutive_regions_overlap_filtered2.sim_snp_indel.bed droSim1_snp_indel_GATC.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa

# sec

#perl get_sequences.pl constitutive_regions_overlap_filtered2.sec_snp_indel.bed droSec1_snpindel_GATC.fa constitutive_regions_overlap_filtered2.sec_snp_indel.fa

##### perform pairwise alignments for each species comparison #####

# zhr and z30

#perl pairwise_aln_FSA.pl zhr z30 constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.z30.fa

# zhr and sim

#perl pairwise_aln_FSA.pl zhr sim constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa

# dm3 and zhr

#perl pairwise_aln_FSA.pl dm3 zhr constitutive_regions_overlap_filtered2.dm3.fa constitutive_regions_overlap_filtered2.zhr.fa

# dm3 and z30

#perl pairwise_aln_FSA.pl dm3 z30 constitutive_regions_overlap_filtered2.dm3.fa constitutive_regions_overlap_filtered2.z30.fa

# z30 and sim

#perl pairwise_aln_FSA.pl z30 sim constitutive_regions_overlap_filtered2.z30.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa

# sim and sec 

#perl pairwise_aln_FSA.pl sim sec constitutive_regions_overlap_filtered2.sim_snp_indel.fa constitutive_regions_overlap_filtered2.sec_snp_indel.fa

##### compare pairwise alignments #####

# zhr and z30

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.zhr.zhr_z30_fsa.fa constitutive_regions_overlap_filtered2.z30.zhr_z30_fsa.fa constitutive_regions_overlap_filtered2.zhr_z30_fsa.SNPs.txt

# dm3 and zhr

perl compare_pairwise.pl constitutive_regions_overlap_filtered2.dm3.dm3_zhr_fsa.fa constitutive_regions_overlap_filtered2.zhr.dm3_zhr_fsa.fa constitutive_regions_overlap_filtered2.dm3_zhr_fsa.SNPs.txt

# dm3 and z30

perl compare_pairwise.pl constitutive_regions_overlap_filtered2.dm3.dm3_z30_fsa.fa constitutive_regions_overlap_filtered2.z30.dm3_z30_fsa.fa constitutive_regions_overlap_filtered2.dm3_z30_fsa.SNPs.txt

# zhr and sim

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.zhr.zhr_sim_fsa.fa constitutive_regions_overlap_filtered2.sim_snp_indel.zhr_sim_fsa.fa constitutive_regions_overlap_filtered2.zhr_sim_fsa.SNPs.txt

# z30 and sim

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.z30.z30_sim_fsa.fa constitutive_regions_overlap_filtered2.sim_snp_indel.z30_sim_fsa.fa constitutive_regions_overlap_filtered2.z30_sim_fsa.SNPs.txt

# sim and sec 

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.sim_snp_indel.sim_sec_fsa.fa constitutive_regions_overlap_filtered2.sec_snp_indel.sim_sec_fsa.fa constitutive_regions_overlap_filtered2.sim_sec_fsa.SNPs.txt