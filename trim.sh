#!/bin/bash
#SBATCH -c8
#SBATCH --mem=32G
#SBATCH --time=8:0:0
#SBATCH -o ./trim

module load trim_galore
module load cutadapt

trim_galore --paired 3n-clone-K.1.fq 3n-clone-K.2.fq --rrbs
