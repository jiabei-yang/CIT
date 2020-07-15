#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcBootIpwSplit.out
#SBATCH -e ../Data/errorRevision/RhcBootIpwSplit.err
#SBATCH --job-name=RhcBootIpwSplit
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD2_RhcBootIpwSplit.R -s 1 -e 1000
