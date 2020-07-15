#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcBootGSplit.out
#SBATCH -e ../Data/errorRevision/RhcBootGSplit.err
#SBATCH --job-name=RhcBootGSplit
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD2_RhcBootGSplit.R -s 1 -e 1000
