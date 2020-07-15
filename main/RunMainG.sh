#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/MainG.out
#SBATCH -e ../Data/errorRevision/MainG.err
#SBATCH --job-name=MainG
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainG.R -s 1 -e 10000
