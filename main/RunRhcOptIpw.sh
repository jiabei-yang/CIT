#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist gpu001
#SBATCH -o ../Data/outputRevision/RhcOptIpw.out
#SBATCH -e ../Data/errorRevision/RhcOptIpw.err
#SBATCH --job-name=IpwOpt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript RhcOptIpw.R 
