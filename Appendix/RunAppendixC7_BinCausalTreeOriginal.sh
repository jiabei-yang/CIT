#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/BinCausalTreeOriginal.out
#SBATCH -e ../Data/errorRevision/BinCausalTreeOriginal.err
#SBATCH --job-name=BinCTOrig
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC7_BinCausalTreeOriginal.R -s 1 -e 1000
