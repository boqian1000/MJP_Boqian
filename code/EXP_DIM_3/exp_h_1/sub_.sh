#!/bin/sh -l
#PBS -l nodes=1:ppn=1
#PBS -l walltime=04:00:00
#PBS -q standby
#PBS -l naccesspolicy=singleuser

cd $PBS_O_WORKDIR
module load python python/2.7.8_intel-16.0.1.150


