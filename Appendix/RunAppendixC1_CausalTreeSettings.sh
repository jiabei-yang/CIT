#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/CausalTreeSettings.out
#SBATCH -e ../Data/errorRevision/CausalTreeSettings.err
#SBATCH --job-name=CTSettings
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC1_CausalTreeSettings.R -s 1 -e 1000
