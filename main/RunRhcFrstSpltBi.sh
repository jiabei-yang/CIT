#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcFrstSpltBi.out
#SBATCH -e ../Data/errorRevision/RhcFrstSpltBi.err
#SBATCH --job-name=RhcFrstSpltBi
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript RhcFrstSpltBi.R 
