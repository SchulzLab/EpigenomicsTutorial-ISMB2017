#!/usr/bin/env sh

##################################################################################################################
# In this script, we will show you how to perform the low level analysis, such as read alignment and peak calling.
# SRR1533847 is a ATAC-seq data for B cell under GEO accession GSE59992
# We assume that all tools have been installed and all necessary data, like referenec genome, have been downloaded. 
##################################################################################################################


#########################################################################################
# We first show the pre-processing of sequence data
#########################################################################################
# downloading of SRA
prefetch -v SRR1533847

# Dump each read into separate file. Files will receive suffix corresponding to read number.
fastq-dump --split-3 ~/ncbi/public/sra/SRR1533847.sra

# Map the reads to reference genome
bowtie2 -x ~/Data/MM10/mm10 -U ./SRR1533847.fastq -S ./B.sam

# Convert sam file to bam file, sort the result and generate the index file
samtools view -Sb ./B.sam > ./B.bam
samtools sort ./B.bam -o ./B_sort.bam
samtools index ./B_sort.bam

# Remove duplicates and bad map quality
samtools rmdup -sS ./B_sort.bam ./B_remdup.bam
samtools view -bq 30 ./B_remdup.bam > ./B.bam
samtools index ./B.bam

# Remove the data unnecessary for downstream analysis
rm B.sam
rm B.fastq
rm B_remdup.bam
rm B_sort.bam
rm B_sort.bam.bai

############################################################################################
# Next, we are going to call the peaks in bam file
############################################################################################

mkdir Peaks
macs2 callpeak -t B.bam -n B --outdir Peaks -f BAM -g mm --nomodel --nolambda --keep-dup all --call-summits

# Merge the overlap peaks and filter bad peak quality
mergeBed -c 9 -o max -i  ./Peaks/B_peaks.narrowPeak > ./Peaks/B_peaks.bed
awk '$4 >= 10' ./Peaks/B_peaks.bed > ./B_peaks.bed

# Remove chrM, and unassembled "random" contigs
sed -i '/chrM/d;/random/d;/chrUn/d' ./B_peaks.bed

##############################################################################################
# Finally, a bam file containing alignment reads and a bed file containing peaks are generated.
# And you can go to footprint calling using these files.
