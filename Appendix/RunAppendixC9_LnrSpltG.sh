#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/LnrSpltG.out
#SBATCH -e ../Data/errorRevision/LnrSpltG.err
#SBATCH --job-name=GLnrSplt
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC9_LnrSpltG.R -s 1 -e 1000 -g 1
