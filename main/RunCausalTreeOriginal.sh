#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute007
#SBATCH -o ../Data/outputRevision/CausalTreeOriginal.out
#SBATCH -e ../Data/errorRevision/CausalTreeOriginal.err
#SBATCH --job-name=CTOrig
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript CausalTreeOriginal.R -s 1 -e 10000
