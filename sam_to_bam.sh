#!/bin/bash
#SBATCH -c1
#SBATCH --mem=32G
#SBATCH --time=2:0:0
#SBATCH -o sam2bam_3n_H

module load samtools

samtools sort 3n_clone_H_meth_02_10.sam > 3n_clone_H_meth_02_10.sorted.bam
samtools index 3n_clone_H_meth_02_10.sorted.bam
