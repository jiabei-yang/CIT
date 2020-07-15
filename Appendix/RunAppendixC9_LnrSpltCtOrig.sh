#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute004
#SBATCH -o ../Data/outputRevision/LnrSpltCtOrig.out
#SBATCH -e ../Data/errorRevision/LnrSpltCtOrig.err
#SBATCH --job-name=CtOrigLnrSplt
#SBATCH --ntasks-per-node=26
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC9_LnrSpltCtOrig.R -s 1 -e 1000 -g 1
