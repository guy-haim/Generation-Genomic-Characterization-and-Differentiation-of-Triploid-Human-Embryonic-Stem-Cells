#!/bin/bash
#SBATCH -c16
#SBATCH --mem=16G
#SBATCH --time=72:0:0
#SBATCH -o ./align_bedgraph

module load bismark
bismark /vol/sci/bio/data/nissim.benvenisty/guyhaim/RRBS -1 3n_clone_L.1_val_1.fq -2 3n_clone_L.2_val_2.fq --multicore 16
bismark_methylation_extractor 3n_clone_L.1_val_1_bismark_bt2.bam --bedGraph
