#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute008
#SBATCH -o ../Data/outputRevision/OutDr.out
#SBATCH -e ../Data/errorRevision/OutDr.err
#SBATCH --job-name=OutDr
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC5_OutDr.R -s 1 -e 1000
