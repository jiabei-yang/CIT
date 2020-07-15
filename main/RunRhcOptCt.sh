#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist gpu001
#SBATCH -o ../Data/outputRevision/RhcOptCt.out
#SBATCH -e ../Data/errorRevision/RhcOptCt.err
#SBATCH --job-name=CtOpt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript RhcOptCt.R 
