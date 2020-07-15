#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute006
#SBATCH -o ../Data/outputRevision/Rf.out
#SBATCH -e ../Data/errorRevision/Rf.err
#SBATCH --job-name=Rf
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC6_Rf.R -s 1 -e 1000
