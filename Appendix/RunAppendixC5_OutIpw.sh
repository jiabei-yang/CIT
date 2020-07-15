#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute001
#SBATCH -o ../Data/outputRevision/OutIpw.out
#SBATCH -e ../Data/errorRevision/OutIpw.err
#SBATCH --job-name=OutIpw
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC5_OutIpw.R -s 1 -e 1000
