#!/usr/bin/env sh

##################################################################################################################
# In this script, we will show you how to perform the low level analysis, such as read alignment and peak calling.
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
bowtie2 -x ~/Data/MM10/mm10 -U ./SRR1533847.fastq -S ./SRR1533847.sam

# Convert sam file to bam file, sort the result and generate the index file
samtools view -Sb ./SRR1533847.sam > ./SRR1533847.bam
samtools sort ./SRR1533847.bam -o ./SRR1533847_sort.bam
samtools index ./SRR1533847_sort.bam

# Remove duplicates and bad map quality
samtools rmdup -sS ./SRR1533847_sort.bam ./SRR1533847_remdup.bam
samtools view -bq 30 ./SRR1533847_remdup.bam > ./SRR1533847.bam
samtools index ./SRR1533847.bam

# Remove the data unnecessary for downstream analysis
rm SRR1533847.sam
rm SRR1533847.fastq
rm SRR1533847_remdup.bam
rm SRR1533847_sort.bam
rm SRR1533847_sort.bam.bai

############################################################################################
# Next, we are going to call the peaks in bam file
############################################################################################

mkdir Peaks
macs2 callpeak -t SRR1533847.bam -n SRR1533847 --outdir Peaks -f BAM -g mm --nomodel --nolambda --keep-dup all --call-summits

# Merge the overlap peaks and filter bad peak quality
mergeBed -c 9 -o max -i  ./Peaks/SRR1533847_peaks.narrowPeak > ./Peaks/SRR1533847_peaks.bed
awk '$4 >= 10' ./Peaks/SRR1533847_peaks.bed > ./SRR1533847_peaks.bed

# Remove chrM, and unassembled "random" contigs
sed -i '/chrM/d;/random/d;/chrUn/d' ./SRR1533847_peaks.bed

##############################################################################################
# Finally, a bam file containing alignment reads and a bed file containing peaks are generated.
# And you can go to footprint calling using these files.
