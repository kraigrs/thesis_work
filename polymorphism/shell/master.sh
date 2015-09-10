#!/bin/sh

##### subtract gaps #####

# zhr and z30

#subtractBed -a constitutive_regions_overlap_filtered2.bed -b Mel_melrsq_gaps_gt19.bed > constitutive_regions_overlap_filtered2_zhr_z30_gapped.bed

# zhr and sim

#subtractBed -a constitutive_regions_overlap_filtered2.bed -b Mel_Simrsq_gaps_gt19.bed > constitutive_regions_overlap_filtered2_zhr_sim_gapped.bed

# sim and sec

#subtractBed -a constitutive_regions_overlap_filtered2.bed -b Sim_sec_gaps_gt19.bed > constitutive_regions_overlap_filtered2_sim_sec_gapped.bed

##### lift exons for each comparison #####

# zhr and z30

#liftOver constitutive_regions_overlap_filtered2_zhr_z30_gapped.bed dm3_zhr_AMBIG_all.chain constitutive_regions_overlap_filtered2_zhr_z30_gapped.zhr.bed constitutive_regions_overlap_filtered2_zhr_z30_gapped.zhr.un

#liftOver constitutive_regions_overlap_filtered2_zhr_z30_gapped.bed dm3_z30_AMBIG_all.chain constitutive_regions_overlap_filtered2_zhr_z30_gapped.z30.bed constitutive_regions_overlap_filtered2_zhr_z30_gapped.z30.un

# zhr and sim

#liftOver constitutive_regions_overlap_filtered2_zhr_sim_gapped.bed dm3_zhr_AMBIG_all.chain constitutive_regions_overlap_filtered2_zhr_sim_gapped.zhr.bed constitutive_regions_overlap_filtered2_zhr_sim_gapped.zhr.un

#liftOver constitutive_regions_overlap_filtered2_zhr_sim_gapped.bed dm3ToDroSim1.over.chain constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim.bed constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim.un

#liftOver constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim.bed droSim1ToDroSim1_snp_indel.over.chain constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim_snp_indel.bed constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim_snp_indel.un

# sim and sec

#liftOver constitutive_regions_overlap_filtered2_sim_sec_gapped.bed dm3ToDroSim1.over.chain constitutive_regions_overlap_filtered2_sim_sec_gapped.sim.bed constitutive_regions_overlap_filtered2_sim_sec_gapped.sim.un

#liftOver constitutive_regions_overlap_filtered2_sim_sec_gapped.sim.bed droSim1ToDroSim1_snp_indel.over.chain constitutive_regions_overlap_filtered2_sim_sec_gapped.sim_snp_indel.bed constitutive_regions_overlap_filtered2_sim_sec_gapped.sim_snp_indel.un

#liftOver constitutive_regions_overlap_filtered2_sim_sec_gapped.bed dm3ToDroSec1.over.chain constitutive_regions_overlap_filtered2_sim_sec_gapped.sec.bed constitutive_regions_overlap_filtered2_sim_sec_gapped.sec.un

#liftOver constitutive_regions_overlap_filtered2_sim_sec_gapped.sec.bed droSec1ToDroSec1_snp_indel.over.chain constitutive_regions_overlap_filtered2_sim_sec_gapped.sec_snp_indel.bed constitutive_regions_overlap_filtered2_sim_sec_gapped.sec_snp_indel.un

##### get lengths of new exons for each comparison #####

# zhr and z30

#perl gene_lengths.pl constitutive_regions_overlap_filtered2_zhr_z30_gapped.bed zhr_z30_exons.txt zhr_z30_genes.txt

# zhr and sim

#perl gene_lengths.pl constitutive_regions_overlap_filtered2_zhr_sim_gapped.bed zhr_sim_exons.txt zhr_sim_genes.txt

# sim and sec

#perl gene_lengths.pl constitutive_regions_overlap_filtered2_sim_sec_gapped.bed sim_sec_exons.txt sim_sec_genes.txt

##### get sequences #####

# zhr and z30

#perl get_sequences.pl constitutive_regions_overlap_filtered2_zhr_z30_gapped.zhr.bed zhr_AMBIG_all.fa constitutive_regions_overlap_filtered2_zhr_z30_gapped.zhr.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2_zhr_z30_gapped.z30.bed z30_AMBIG_all.fa constitutive_regions_overlap_filtered2_zhr_z30_gapped.z30.fa

# zhr and sim

#perl get_sequences.pl constitutive_regions_overlap_filtered2_zhr_sim_gapped.zhr.bed zhr_AMBIG_all.fa constitutive_regions_overlap_filtered2_zhr_sim_gapped.zhr.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim_snp_indel.bed droSim1_snp_indel_GATC.fa constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim_snp_indel.fa

# sim and sec

#perl get_sequences.pl constitutive_regions_overlap_filtered2_sim_sec_gapped.sim_snp_indel.bed droSim1_snp_indel_GATC.fa constitutive_regions_overlap_filtered2_sim_sec_gapped.sim_snp_indel.fa

#perl get_sequences.pl constitutive_regions_overlap_filtered2_sim_sec_gapped.sec_snp_indel.bed droSec1_snpindel_GATC.fa constitutive_regions_overlap_filtered2_sim_sec_gapped.sec_snp_indel.fa

##### perform pairwise alignments for each species comparison #####

# zhr and z30

#perl pairwise_aln_FSA.pl constitutive_regions_overlap_filtered2_zhr_z30_gapped.zhr.fa constitutive_regions_overlap_filtered2_zhr_z30_gapped.z30.fa

# zhr and sim

#perl pairwise_aln_FSA.pl constitutive_regions_overlap_filtered2_zhr_sim_gapped.zhr.fa constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim_snp_indel.fa

# sim and sec 

#perl pairwise_aln_FSA.pl constitutive_regions_overlap_filtered2_sim_sec_gapped.sim_snp_indel.fa constitutive_regions_overlap_filtered2_sim_sec_gapped.sec_snp_indel.fa

##### compare pairwise alignments #####

# zhr and z30

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2_zhr_z30_gapped.zhr_fsa.fa constitutive_regions_overlap_filtered2_zhr_z30_gapped.z30_fsa.fa constitutive_regions_overlap_filtered2_zhr_z30_gapped.SNPs.txt

# zhr and sim

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2_zhr_sim_gapped.zhr_fsa.fa constitutive_regions_overlap_filtered2_zhr_sim_gapped.sim_snp_indel_fsa.fa constitutive_regions_overlap_filtered2_zhr_sim_gapped.SNPs.txt

# sim and sec 

#perl compare_pairwise.pl constitutive_regions_overlap_filtered2_sim_sec_gapped.sim_snp_indel_fsa.fa constitutive_regions_overlap_filtered2_sim_sec_gapped.sec_snp_indel_fsa.fa constitutive_regions_overlap_filtered2_sim_sec_gapped.SNPs.txt