#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist gpu001
#SBATCH -o ../Data/outputRevision/RhcOptG.out
#SBATCH -e ../Data/errorRevision/RhcOptG.err
#SBATCH --job-name=GOpt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript RhcOptG.R 
