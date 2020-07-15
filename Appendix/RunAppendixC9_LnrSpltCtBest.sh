#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/LnrSpltCtBest.out
#SBATCH -e ../Data/errorRevision/LnrSpltCtBest.err
#SBATCH --job-name=CtBestLnrSplt
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC9_LnrSpltCtBest.R -s 1 -e 1000 -g 1
