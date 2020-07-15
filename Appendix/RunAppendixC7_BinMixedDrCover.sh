#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/BinMixedDrCover.out
#SBATCH -e ../Data/errorRevision/BinMixedDrCover.err
#SBATCH --job-name=BinDrCover
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC7_BinMixedDrCover.R -s 1 -e 1000
