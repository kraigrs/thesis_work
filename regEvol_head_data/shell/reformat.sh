#!/bin/sh

# mel_mel

FILES=../../mel_mel_data/*.sam
for f in $FILES
do
  var1=$f qsub samToBed.sh
done

# sim_sec

FILES=../../sim_sec_data/*.sam
for f in $FILES
do
  var1=$f qsub samToBed.sh
done

# mel_sim

FILES=../../mel_sim_data/*.sam
for f in $FILES
do
  var1=$f qsub samToBed.sh
done

# mel_sec

FILES=../../mel_sec_data/*.sam
for f in $FILES
do
  var1=$f qsub samToBed.sh
done