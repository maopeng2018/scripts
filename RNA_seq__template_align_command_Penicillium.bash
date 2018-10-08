#!/bin/bash
bwa index -p ../Genomes_assembled/Penicillium-assembly_bwa_ref ../Genomes_assembled/Penicillium-assembly.fa
/Data/genome_assembly/RSEM-master/rsem-prepare-reference --gtf ../Genomes_assembled/annotation-Penicillium.gtf --bowtie2 ../Genomes_assembled/Penicillium-assembly.fa ../Genomes_assembled/Penicillium_bowtie2_ref


bwa aln ../Genomes_assembled/Penicillium-assembly_bwa_ref  7047-03-09_S1_L001_R1_001.fastq.gz >7047-03-09_S1_L001_R1_001_bwa.sai -t 18
bwa aln ../Genomes_assembled/Penicillium-assembly_bwa_ref  7047-03-09_S1_L001_R2_001.fastq.gz >7047-03-09_S1_L001_R2_001_bwa.sai -t 18
bwa sampe ../Genomes_assembled/Penicillium-assembly_bwa_ref 7047-03-09_S1_L001_R1_001_bwa.sai 7047-03-09_S1_L001_R2_001_bwa.sai 7047-03-09_S1_L001_R1_001.fastq.gz 7047-03-09_S1_L001_R2_001.fastq.gz >7047-03-09_S1_L001_001_bwa.sam 
samtools view -bS 7047-03-09_S1_L001_001_bwa.sam >7047-03-09_S1_L001_001_bwa.bam
samtools sort 7047-03-09_S1_L001_001_bwa.bam 7047-03-09_S1_L001_001_bwa_sorted
samtools index 7047-03-09_S1_L001_001_bwa_sorted.bam
samtools flagstat 7047-03-09_S1_L001_001_bwa.bam >7047-03-09_S1_L001_001_mapped_bwa_summary

bowtie2 -p 18 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --no-mixed --no-discordant -x ../Genomes_assembled/Penicillium_bowtie2_ref -1 7047-03-09_S1_L001_R1_001.fastq.gz -2 7047-03-09_S1_L001_R2_001.fastq.gz -S 7047-03-09_S1_L001_001_bowtie2.sam                          
samtools view -bS 7047-03-09_S1_L001_001_bowtie2.sam >7047-03-09_S1_L001_001_bowtie2.bam 
samtools sort 7047-03-09_S1_L001_001_bowtie2.bam 7047-03-09_S1_L001_001_bowtie2_sorted
samtools index 7047-03-09_S1_L001_001_bowtie2_sorted.bam 
samtools flagstat 7047-03-09_S1_L001_001_bowtie2.bam >7047-03-09_S1_L001_001_mapped_bowtie2_summary 
 
/Data/genome_assembly/RSEM-master/rsem-calculate-expression --no-bam-output -p 18 --bam --paired-end ./7047-03-09_S1_L001_001_bowtie2.bam ../Genomes_assembled/Penicillium_bowtie2_ref 7047-03-09_S1_L001_001_bowtie2F                  



