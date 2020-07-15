#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/BinMixedDr.out
#SBATCH -e ../Data/errorRevision/BinMixedDr.err
#SBATCH --job-name=BinDr
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC7_BinMixedDr.R -s 1 -e 1000
