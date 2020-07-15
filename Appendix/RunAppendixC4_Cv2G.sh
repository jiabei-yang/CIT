#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute006
#SBATCH -o ../Data/outputRevision/Cv2G.out
#SBATCH -e ../Data/errorRevision/Cv2G.err
#SBATCH --job-name=Cv2G
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC4_Cv2G.R -s 1 -e 1000
