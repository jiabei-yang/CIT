#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute003
#SBATCH -o ../Data/outputRevision/SplitDr.out
#SBATCH -e ../Data/errorRevision/SplitDr.err
#SBATCH --job-name=SplitDr
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC5_SplitDr.R -s 1 -e 1000
