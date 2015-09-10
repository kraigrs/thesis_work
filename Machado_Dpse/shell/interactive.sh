#!/bin/sh

qsub -I -V -A lsa_flux -l qos=flux,nodes=1:ppn=1,pmem=4000mb,walltime=1:00:00 -q flux