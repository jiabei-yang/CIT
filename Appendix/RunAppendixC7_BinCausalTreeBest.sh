#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/BinCausalTreeBest.out
#SBATCH -e ../Data/errorRevision/BinCausalTreeBest.err
#SBATCH --job-name=BinCTBest
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC7_BinCausalTreeBest.R -s 1 -e 1000
