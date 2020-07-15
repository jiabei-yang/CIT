#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcBootDrSplit.out
#SBATCH -e ../Data/errorRevision/RhcBootDrSplit.err
#SBATCH --job-name=RhcBootDrSplit
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD2_RhcBootDrSplit.R -s 1 -e 1000
