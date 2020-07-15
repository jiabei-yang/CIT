#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist gpu001
#SBATCH -o ../Data/outputRevision/LnrSpltCtSettings.out
#SBATCH -e ../Data/errorRevision/LnrSpltCtSettings.err
#SBATCH --job-name=CtSetsLnrSplt
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC9_LnrSpltCtSettings.R -s 1 -e 1000 -g 1
