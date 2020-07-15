#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute002
#SBATCH -o ../Data/outputRevision/BinCausalTreeSettings.out
#SBATCH -e ../Data/errorRevision/BinCausalTreeSettings.err
#SBATCH --job-name=BinCTSettings
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC7_BinCausalTreeSettings.R -s 1 -e 1000
