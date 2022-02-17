#!/bin/bash
#SBATCH -c32
#SBATCH --mem=64G
#SBATCH --time=72:0:0
#SBATCH -o methylation_3n_K

module load methyldackel
MethylDackel extract GRCh38.primary_assembly.genome.fa /vol/sci/bio/data/nissim.benvenisty/guyhaim/bwameth_env/3n_clone_K_meth_02_10.sorted.bam
