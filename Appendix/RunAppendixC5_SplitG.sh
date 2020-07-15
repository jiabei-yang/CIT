#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute006
#SBATCH -o ../Data/outputRevision/SplitG.out
#SBATCH -e ../Data/errorRevision/SplitG.err
#SBATCH --job-name=SplitG
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC5_SplitG.R -s 1 -e 1000
