#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute002
#SBATCH -o ../Data/outputRevision/SmSigAr1CtSettings.out
#SBATCH -e ../Data/errorRevision/SmSigAr1CtSettings.err
#SBATCH --job-name=SmCTSettings
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC8_SmSigAr1CtSettings.R -s 1 -e 1000 -g 1
