#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute002
#SBATCH -o ../Data/outputRevision/Cv2Ipw.out
#SBATCH -e ../Data/errorRevision/Cv2Ipw.err
#SBATCH --job-name=Cv2Ipw
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript RunAppendixC4_Cv2Ipw.R -s 1 -e 1000
