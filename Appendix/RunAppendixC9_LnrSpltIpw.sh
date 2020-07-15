#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute002
#SBATCH -o ../Data/outputRevision/LnrSpltIpw.out
#SBATCH -e ../Data/errorRevision/LnrSpltIpw.err
#SBATCH --job-name=IpwLnrSplt
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC9_LnrSpltIpw.R -s 1 -e 1000 -g 1
