#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcSpltRootCi.out
#SBATCH -e ../Data/errorRevision/RhcSpltRootCi.err
#SBATCH --job-name=RhcSpltRootCi
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD3_RhcSpltRootCi.R -s 1 -e 1000
