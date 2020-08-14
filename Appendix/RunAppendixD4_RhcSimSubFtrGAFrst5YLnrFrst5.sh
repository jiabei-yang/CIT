#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/outputRevision/RhcSimSubFtrGAFrst5YLnrFrst5.out
#SBATCH -e ../Data/errorRevision/RhcSimSubFtrGAFrst5YLnrFrst5.err
#SBATCH --job-name=RhcSimSubFtrGAFrst5YLnrFrst5
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixD4_RhcSimSubFtrGAFrst5YLnrFrst5.R -s 1 -e 1000
