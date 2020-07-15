#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute006
#SBATCH -o ../Data/outputRevision/LnrSpltDr.out
#SBATCH -e ../Data/errorRevision/LnrSpltDr.err
#SBATCH --job-name=DrLnrSplt
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC9_LnrSpltDr.R -s 1 -e 1000 -g 1
