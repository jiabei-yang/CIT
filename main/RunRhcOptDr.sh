#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist gpu001
#SBATCH -o ../Data/outputRevision/RhcOptDr.out
#SBATCH -e ../Data/errorRevision/RhcOptDr.err
#SBATCH --job-name=DrOpt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript RhcOptDr.R 
