#!/bin/sh

# lift introns and regulatory regions
 
liftOver ../../references/dm3_regulatory_regions_merged.bed ../../genomes/zhr/dm3_zhr_AMBIG_all.chain ../../references/zhr_z30/dm3_regulatory_regions_merged.zhr.bed ../../references/zhr_z30/dm3_regulatory_regions_merged.zhr.un -minMatch=0.5

liftOver ../../references/dm3_regulatory_regions_merged.bed ../../genomes/z30/dm3_z30_AMBIG_all.chain ../../references/zhr_z30/dm3_regulatory_regions_merged.z30.bed ../../references/zhr_z30/dm3_regulatory_regions_merged.z30.un -minMatch=0.5

liftOver ../../references/dm3_regulatory_regions_merged.bed ../../genomes/dm3/dm3ToDroSim1.over.chain ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1.bed ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1.un -minMatch=0.5

liftOver ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1.bed ../../genomes/tsimbazaza/droSim1ToDroSim1_snp_indel.over.chain ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1_snp_indel.bed ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1_snp_indel.un -minMatch=0.5

liftOver ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1.un ../../genomes/dm3/dm3ToSim_extra.over.chain ../../references/zhr_tsim/dm3_regulatory_regions_merged.Sim_extra.bed ../../references/zhr_tsim/dm3_regulatory_regions_merged.Sim_extra.un -minMatch=0.5

liftOver ../../references/dm3_regulatory_regions_merged.bed ../../genomes/dm3/dm3ToDroSec1.over.chain ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.bed ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.un -minMatch=0.5

liftOver ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.bed ../../genomes/droSec1/droSec1ToDroSec1_snp_indel.over.chain ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1_snp_indel.bed ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1_snp_indel.un -minMatch=0.5

liftOver ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.un ../../genomes/dm3/dm3ToSec_extra.over.chain ../../references/tsim_droSec1/dm3_regulatory_regions_merged.Sec_extra.bed ../../references/tsim_droSec1/dm3_regulatory_regions_merged.Sec_extra.un -minMatch=0.5

# get sequences

perl ../perl/get_sequences.pl ../../references/zhr_z30/dm3_regulatory_regions_merged.zhr.bed ../../genomes/zhr/zhr_AMBIG_all.fa ../../references/zhr_z30/dm3_regulatory_regions_merged.zhr.fa

perl ../perl/get_sequences.pl ../../references/zhr_z30/dm3_regulatory_regions_merged.z30.bed ../../genomes/z30/z30_AMBIG_all.fa ../../references/zhr_z30/dm3_regulatory_regions_merged.z30.fa

perl ../perl/get_sequences.pl ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1_snp_indel.bed ../../genomes/tsimbazaza/droSim1_snp_indel.fa ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1_snp_indel.fa

perl ../perl/get_sequences.pl ../../references/zhr_tsim/dm3_regulatory_regions_merged.Sim_extra.bed ../../genomes/tsimbazaza/Sim_k35_extra.fa ../../references/zhr_tsim/dm3_regulatory_regions_merged.Sim_extra.fa

cat ../../references/zhr_tsim/dm3_regulatory_regions_merged.droSim1_snp_indel.fa ../../references/zhr_tsim/dm3_regulatory_regions_merged.Sim_extra.fa > ../../references/zhr_tsim/dm3_regulatory_regions_merged.tsim.fa

perl ../perl/get_sequences.pl ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1_snp_indel.bed ../../genomes/droSec1/droSec1_snp_indel.fa ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1_snp_indel.fa

perl ../perl/get_sequences.pl ../../references/tsim_droSec1/dm3_regulatory_regions_merged.Sec_extra.bed ../../genomes/droSec1/Sec_rsq_k35_extra.fa ../../references/tsim_droSec1/dm3_regulatory_regions_merged.Sec_extra.fa

cat ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1_snp_indel.fa ../../references/tsim_droSec1/dm3_regulatory_regions_merged.Sec_extra.fa > ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.fa

# perform pairwise alignment

perl ../perl/pairwise_aln_FSA.pl zhr z30 ../../references/zhr_z30/dm3_regulatory_regions_merged.zhr.fa ../../references/zhr_z30/dm3_regulatory_regions_merged.z30.fa

perl ../perl/pairwise_aln_FSA.pl tsim droSec1 ../../references/tsim_droSec1/dm3_regulatory_regions_merged.tsim.fa ../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.fa

perl ../perl/pairwise_aln_FSA.pl zhr tsim ../../references/zhr_tsim/dm3_regulatory_regions_merged.zhr.fa ../../references/zhr_tsim/dm3_regulatory_regions_merged.tsim.fa

# compare alignments to find differentiating sites

#perl ../perl/compare_pairwise.pl ../../validation/zhr_candidate_regs.zhr_z30_fsa.fa ../../validation/zhr_candidate_regs.zhr_z30_fsa.fa ../../validation/zhr_z30_fsa.SNPs.txt