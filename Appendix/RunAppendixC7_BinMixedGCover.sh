#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute003
#SBATCH -o ../Data/outputRevision/BinMixedGCover.out
#SBATCH -e ../Data/errorRevision/BinMixedGCover.err
#SBATCH --job-name=BinGCover
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC7_BinMixedGCover.R -s 1 -e 1000
