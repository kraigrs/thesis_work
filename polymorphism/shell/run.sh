#!/bin/sh


##### lift exons for each comparison #####

# zhr

liftOver constitutive_regions_overlap_filtered2.bed dm3_zhr_AMBIG_all.chain constitutive_regions_overlap_filtered2.zhr.bed constitutive_regions_overlap_filtered2.zhr.un 

# z30

liftOver constitutive_regions_overlap_filtered2.bed dm3_z30_AMBIG_all.chain constitutive_regions_overlap_filtered2.z30.bed constitutive_regions_overlap_filtered2.z30.un 

# sim

liftOver constitutive_regions_overlap_filtered2.bed dm3ToDroSim1.over.chain constitutive_regions_overlap_filtered2.sim.bed constitutive_regions_overlap_filtered2.sim.un 

liftOver constitutive_regions_overlap_filtered2.sim.bed droSim1ToDroSim1_snp_indel.over.chain constitutive_regions_overlap_filtered2.sim_snp_indel.bed constitutive_regions_overlap_filtered2.sim_snp_indel.un 

liftOver constitutive_regions_overlap_filtered2.sim.un dm3ToSim_extra.over.chain constitutive_regions_overlap_filtered2.sim_extra.bed constitutive_regions_overlap_filtered2.sim_extra.un 

# sec

liftOver constitutive_regions_overlap_filtered2.bed dm3ToDroSec1.over.chain constitutive_regions_overlap_filtered2.sec.bed constitutive_regions_overlap_filtered2.sec.un 

liftOver constitutive_regions_overlap_filtered2.sec.bed droSec1ToDroSec1_snp_indel.over.chain constitutive_regions_overlap_filtered2.sec_snp_indel.bed constitutive_regions_overlap_filtered2.sec_snp_indel.un 

liftOver constitutive_regions_overlap_filtered2.sec.un dm3ToSec_extra.over.chain constitutive_regions_overlap_filtered2.sec_extra.bed constitutive_regions_overlap_filtered2.sec_extra.un 


##### get sequences #####

# zhr

perl get_sequences.pl constitutive_regions_overlap_filtered2.zhr.bed zhr_AMBIG_all.fa constitutive_regions_overlap_filtered2.zhr.fa

# z30

perl get_sequences.pl constitutive_regions_overlap_filtered2.z30.bed z30_AMBIG_all.fa constitutive_regions_overlap_filtered2.z30.fa

# sim

perl get_sequences.pl constitutive_regions_overlap_filtered2.sim_snp_indel.bed droSim1_snp_indel_GATC.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa

perl get_sequences.pl constitutive_regions_overlap_filtered2.sim_extra.bed Sim_k35_extra.fa constitutive_regions_overlap_filtered2.sim_extra.fa

cat constitutive_regions_overlap_filtered2.sim_snp_indel.fa constitutive_regions_overlap_filtered2.sim_extra.fa > constitutive_regions_overlap_filtered2.sim_snp_indel_extra.fa

# sec

perl get_sequences.pl constitutive_regions_overlap_filtered2.sec_snp_indel.bed droSec1_snpindel_GATC.fa constitutive_regions_overlap_filtered2.sec_snp_indel.fa

perl get_sequences.pl constitutive_regions_overlap_filtered2.sec_extra.bed Sec_rsq_k35_extra.fa constitutive_regions_overlap_filtered2.sec_extra.fa

cat constitutive_regions_overlap_filtered2.sec_snp_indel.fa constitutive_regions_overlap_filtered2.sec_extra.fa > constitutive_regions_overlap_filtered2.sec_snp_indel_extra.fa


##### perform pairwise alignments for each species comparison #####

# zhr and z30

#perl pairwise_aln_FSA.pl zhr z30 constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.z30.fa

# zhr and sim

#perl pairwise_aln_FSA.pl zhr sim constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.sim_snp_indel_extra.fa

# sim and sec 

#perl pairwise_aln_FSA.pl sim sec constitutive_regions_overlap_filtered2.sim_snp_indel_extra.fa constitutive_regions_overlap_filtered2.sec_snp_indel_extra.fa
