#!/bin/sh

qsub -I -V -A lsa_flux -l procs=1,qos=flux,walltime=1:00:00 -q flux