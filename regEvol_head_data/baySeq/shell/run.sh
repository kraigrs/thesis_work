#!/bin/sh

var1=./mel_mel/baySeq_mel_mel1_female qsub baySeq.sh
var1=./mel_mel/baySeq_mel_mel1_male qsub baySeq.sh
var1=./mel_mel/baySeq_mel_mel2_female qsub baySeq.sh
var1=./mel_mel/baySeq_mel_mel2_male qsub baySeq.sh
var1=./mel_mel/baySeq_mel_mel3_female qsub baySeq.sh
var1=./mel_mel/baySeq_mel_mel3_male qsub baySeq.sh

var1=./sim_sec/baySeq_sim_sec_female qsub baySeq.sh
var1=./sim_sec/baySeq_sim_sec_male qsub baySeq.sh

var1=./mel_sim/baySeq_mel_sim_female qsub baySeq.sh
var1=./mel_sim/baySeq_mel_sim_male qsub baySeq.sh

var1=./mel_sec/baySeq_mel_sec_female qsub baySeq.sh
