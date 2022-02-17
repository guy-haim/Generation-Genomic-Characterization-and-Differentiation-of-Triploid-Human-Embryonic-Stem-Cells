#!/bin/bash

module load elkind


#zcat /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/all_FASTQ_files/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_L001_R1_001.fastq.gz# /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/all_FASTQ_files/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_L002_R1_001.fastq.gz# /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/all_FASTQ_files/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_L003_R1_001.fastq.gz# /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/all_FASTQ_files/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_L004_R1_001.fastq.gz > /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/fastq/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H.fastq

trimmomatic SE /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/fastq/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H.fastq /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/trim/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_trimmed.fastq ILLUMINACLIP:/vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#!/bin/bash
module load r4
module load samtools/1.10
module load bio

# Human:
STAR --runThreadN 8 --genomeDir /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/genome_index --readFilesIn /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/trim/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_trimmed.fastq --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattributes NM --outFileNamePrefix /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/star_out/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_

# Mouse:
STAR --runThreadN 8 --genomeDir /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/genome_index_mouse --readFilesIn /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/trim/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_trimmed.fastq --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattributes NM --outFileNamePrefix /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/star_mouse_out/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_mouse_

human=/vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/star_out/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_Aligned.sortedByCoord.out.bam
mouse=/vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/star_mouse_out/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_mouse_Aligned.sortedByCoord.out.bam
destination=/vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/bam/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H

Rscript --vanilla ./xenofilterr.R ${human} ${mouse} ${destination}

featureCounts -T 8 -a /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/gencode.v35.primary_assembly.annotation.gtf -o /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/feature_counts/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_counts.txt /vol/sci/bio/data/nissim.benvenisty/guyhaim/Xeno_transcriptome_alignment/bam/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H/Filtered_bams/T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_Aligned.sortedByCoord.out_Filtered.bam

