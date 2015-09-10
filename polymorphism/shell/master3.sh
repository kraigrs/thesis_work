#!/bin/sh

##### lift exons for each comparison #####

# zhr and z30

#liftOver constitutive_regions_overlap_filtered2.zhr.un dm3_zhr_AMBIG_all.chain constitutive_regions_overlap_filtered2.zhr_failed.bed constitutive_regions_overlap_filtered2.zhr_failed.un -minMatch=0.5

#liftOver constitutive_regions_overlap_filtered2.z30.un dm3_z30_AMBIG_all.chain constitutive_regions_overlap_filtered2.z30_failed.bed constitutive_regions_overlap_filtered2.z30_failed.un -minMatch=0.5

# sim

#liftOver constitutive_regions_overlap_filtered2.sim.un dm3ToDroSim1.over.chain constitutive_regions_overlap_filtered2.sim_failed.bed constitutive_regions_overlap_filtered2.sim_failed.un -minMatch=0.5

#liftOver constitutive_regions_overlap_filtered2.sim_failed.bed droSim1ToDroSim1_snp_indel.over.chain constitutive_regions_overlap_filtered2.sim_snp_indel_failed.bed constitutive_regions_overlap_filtered2.sim_snp_indel_failed.un -minMatch=0.5

#liftOver constitutive_regions_overlap_filtered2.sim.un dm3ToSim_extra.over.chain constitutive_regions_overlap_filtered2.sim_extra.bed constitutive_regions_overlap_filtered2.sim_extra.un -minMatch=0.5

#liftOver constitutive_regions_overlap_filtered2.sim_failed.un dm3ToSim_extra.over.chain constitutive_regions_overlap_filtered2.sim_failed_extra.bed constitutive_regions_overlap_filtered2.sim_failed_extra.un -minMatch=0.5

# sec

liftOver constitutive_regions_overlap_filtered2.sec.un dm3ToDroSec1.over.chain constitutive_regions_overlap_filtered2.sec_failed.bed constitutive_regions_overlap_filtered2.sec_failed.un -minMatch=0.5

liftOver constitutive_regions_overlap_filtered2.sec_failed.bed droSec1ToDroSec1_snp_indel.over.chain constitutive_regions_overlap_filtered2.sec_snp_indel_failed.bed constitutive_regions_overlap_filtered2.sec_snp_indel_failed.un -minMatch=0.5

liftOver constitutive_regions_overlap_filtered2.sec.un dm3ToSec_extra.over.chain constitutive_regions_overlap_filtered2.sec_extra.bed constitutive_regions_overlap_filtered2.sec_extra.un -minMatch=0.5

##### get sequences #####

# dm3

#perl get_sequences.pl constitutive_regions_overlap_filtered2.bed ../Dmel/dm3.fa constitutive_regions_overlap_filtered2.dm3.fa

# zhr and z30

#perl get_sequences.pl constitutive_regions_overlap_filtered2.zhr.bed zhr_AMBIG_all.fa constitutive_regions_overlap_filtered2.zhr.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2.z30_failed.bed z30_AMBIG_all.fa constitutive_regions_overlap_filtered2.z30_failed.fa

# sim

#perl get_sequences.pl constitutive_regions_overlap_filtered2.sim_snp_indel.bed droSim1_snp_indel_GATC.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2.sim_extra.bed Sim_k35_extra.fa constitutive_regions_overlap_filtered2.sim_extra.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2.sim_failed_extra.bed Sim_k35_extra.fa constitutive_regions_overlap_filtered2.sim_failed_extra.fa

# sec

#perl get_sequences.pl constitutive_regions_overlap_filtered2.sec_snp_indel.bed droSec1_snpindel_GATC.fa constitutive_regions_overlap_filtered2.sec_snp_indel.fa

perl get_sequences.pl constitutive_regions_overlap_filtered2.sec_extra.bed Sec_rsq_k35_extra.fa constitutive_regions_overlap_filtered2.sec_extra.fa

##### perform pairwise alignments for each species comparison #####

# zhr and z30

#perl pairwise_aln_FSA.pl zhr z30 constitutive_regions_overlap_filtered2.zhr_failed.fa constitutive_regions_overlap_filtered2.z30_failed.fa

# zhr and sim

#perl pairwise_aln_FSA.pl zhr sim constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa

perl pairwise_aln_FSA.pl zhr sim_failed_extra constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.sim_failed_extra.fa

# sim and sec 

#perl pairwise_aln_FSA.pl sim sec constitutive_regions_overlap_filtered2.sim_snp_indel.fa constitutive_regions_overlap_filtered2.sec_snp_indel.fa

#perl pairwise_aln_FSA.pl sim sec constitutive_regions_overlap_filtered2.sim_snp_indel_failed.fa constitutive_regions_overlap_filtered2.sec_snp_indel_failed.fa

#perl pairwise_aln_FSA.pl sim sec constitutive_regions_overlap_filtered2.sim_snp_indel_failed.fa constitutive_regions_overlap_filtered2.sec_snp_indel_failed.fa

perl pairwise_aln_FSA.pl sim sec constitutive_regions_overlap_filtered2.sim_snp_indel_extra.fa constitutive_regions_overlap_filtered2.sec_snp_indel_extra.fa

##### compare pairwise alignments #####

# zhr and z30

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.zhr.zhr_z30_fsa.fa constitutive_regions_overlap_filtered2.z30.zhr_z30_fsa.fa constitutive_regions_overlap_filtered2.zhr_z30_fsa.SNPs.txt

# zhr and sim

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.zhr.zhr_sim_fsa.fa constitutive_regions_overlap_filtered2.sim_snp_indel.zhr_sim_fsa.fa constitutive_regions_overlap_filtered2.zhr_sim_fsa.SNPs.txt

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.zhr.zhr_sim_extra_fsa.fa constitutive_regions_overlap_filtered2.sim_extra.zhr_sim_extra_fsa.fa

perl compare_pairwise.pl constitutive_regions_overlap_filtered2.zhr.zhr_sim_failed_extra_fsa.fa constitutive_regions_overlap_filtered2.sim_failed_extra.zhr_sim_failed_extra_fsa.fa

# sim and sec 

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2.sim_snp_indel.sim_sec_fsa.fa constitutive_regions_overlap_filtered2.sec_snp_indel.sim_sec_fsa.fa constitutive_regions_overlap_filtered2.sim_sec_fsa.SNPs.txt

perl compare_pairwise.pl constitutive_regions_overlap_filtered2.sim_snp_indel_extra.sim_sec_fsa.fa constitutive_regions_overlap_filtered2.sec_snp_indel_extra.sim_sec_fsa.fa