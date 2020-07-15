#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcSimSubFtrDrAFrst5YLnrFrst5.out
#SBATCH -e ../Data/errorRevision/RhcSimSubFtrDrAFrst5YLnrFrst5.err
#SBATCH --job-name=RhcSimSubFtrDrAFrst5YLnrFrst5
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD4_RhcSimSubFtrDrAFrst5YLnrFrst5.R -s 1 -e 1000
