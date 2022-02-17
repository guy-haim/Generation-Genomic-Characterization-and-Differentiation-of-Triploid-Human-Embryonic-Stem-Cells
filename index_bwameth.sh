#!/bin/bash
#SBATCH -c16
#SBATCH --mem=32G
#SBATCH --time=12:0:0
#SBATCH -o output_index

. /vol/sci/bio/data/nissim.benvenisty/guyhaim/bwameth_env/bin/activate

bwameth.py index /vol/sci/bio/data/nissim.benvenisty/guyhaim/RRBS/GRCh38.primary_assembly.genome.fa
