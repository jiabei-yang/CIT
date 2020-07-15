#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute003
#SBATCH -o ../Data/outputRevision/SplitIpw.out
#SBATCH -e ../Data/errorRevision/SplitIpw.err
#SBATCH --job-name=SplitIpw
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC5_SplitIpw.R -s 1 -e 1000
