#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/RhcSpltFrstSpltCi.out
#SBATCH -e ../Data/errorRevision/RhcSpltFrstSpltCi.err
#SBATCH --job-name=RhcFrstSpltCi
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD3_RhcSpltFrstSpltCi.R 
