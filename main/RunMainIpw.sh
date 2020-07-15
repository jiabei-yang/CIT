#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute007
#SBATCH -o ../Data/outputRevision/MainIpw.out
#SBATCH -e ../Data/errorRevision/MainIpw.err
#SBATCH --job-name=MainIpw
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainIpw.R -s 1 -e 10000
