#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute008
#SBATCH -o ../Data/outputRevision/CausalTreeBest.out
#SBATCH -e ../Data/errorRevision/CausalTreeBest.err
#SBATCH --job-name=CTBest
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript CausalTreeBest.R -s 1 -e 10000
