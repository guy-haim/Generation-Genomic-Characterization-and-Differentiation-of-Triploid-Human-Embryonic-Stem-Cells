#!/bin/bash
#SBATCH -c32
#SBATCH --mem=64G
#SBATCH --time=90:0:0
#SBATCH -o methylation_1n


. /vol/sci/bio/data/nissim.benvenisty/guyhaim/bwameth_env/bin/activate

bwameth.py --threads 32 \
     --reference /vol/sci/bio/data/nissim.benvenisty/guyhaim/RRBS/GRCh38.primary_assembly.genome.fa \
     /vol/sci/bio/data/nissim.benvenisty/guyhaim/RRBS/3n-clone-H.1_val_1.fq /vol/sci/bio/data/nissim.benvenisty/guyhaim/RRBS/3n-clone-H.2_val_2.fq > 3n_clone_H_meth_02_10.sam
