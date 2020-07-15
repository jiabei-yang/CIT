#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute005
#SBATCH -o ../Data/outputRevision/MainDr.out
#SBATCH -e ../Data/errorRevision/MainDr.err
#SBATCH --job-name=MainDr
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainDr.R -s 1 -e 10000
