#!/bin/sh

var1=../../genomes/lifted/zhr_AMBIG_all.bowtie_v0_m1.bed var2=../../genomes/lifted/constitutive_regions_overlap_filtered2.zhr.bed var3=zhr qsub coverageBed.sh

var1=../../genomes/lifted/z30_AMBIG_all.bowtie_v0_m1.bed var2=../../genomes/lifted/constitutive_regions_overlap_filtered2.z30.bed var3=z30 qsub coverageBed.sh

var1=../../genomes/lifted/droSim1_snp_indel.bowtie_v0_m1.bed var2=../../genomes/lifted/constitutive_regions_overlap_filtered2.sim.bed var3=sim qsub coverageBed.sh

var1=../../genomes/lifted/droSec1_snp_indel.bowtie_v0_m1.bed var2=../../genomes/lifted/constitutive_regions_overlap_filtered2.sec.bed var3=sec qsub coverageBed.sh

var1=../../genomes/lifted/sim_extra.bowtie_v0_m1.bed var2=../../genomes/lifted/constitutive_regions_overlap_filtered2.sim_extra.bed var3=sim_extra qsub coverageBed.sh

var1=../../genomes/lifted/sec_extra.bowtie_v0_m1.bed var2=../../genomes/lifted/constitutive_regions_overlap_filtered2.sec_extra.bed var3=sec_extra qsub coverageBed.sh
