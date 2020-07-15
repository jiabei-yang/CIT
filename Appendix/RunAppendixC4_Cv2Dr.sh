#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodelist compute009
#SBATCH -o ../Data/outputRevision/Cv2Dr.out
#SBATCH -e ../Data/errorRevision/Cv2Dr.err
#SBATCH --job-name=Cv2Dr
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript AppendixC4_Cv2Dr.R -s 1 -e 1000
