#!/bin/bash
#SBATCH -c8
#SBATCH --mem=32G
#SBATCH --time=5:0:0
#SBATCH -o ./bismark

module load bismark
bismark_genome_preparation /vol/sci/bio/data/nissim.benvenisty/guyhaim/RRBS
