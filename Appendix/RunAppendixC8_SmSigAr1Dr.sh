#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute006
#SBATCH -o ../Data/outputRevision/SmSigAr1Dr.out
#SBATCH -e ../Data/errorRevision/SmSigAr1Dr.err
#SBATCH --job-name=DrSmSigAr1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC8_SmSigAr1Dr.R -s 1 -e 1000 -g 1
