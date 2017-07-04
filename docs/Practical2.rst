==================================================================
Practical II - Footprint calling & Transcription factor prediction
==================================================================

Step1: Footprint calling
-----------------------------------------------

After installation of HINT, go to the *EpigenomicsTutorial-ISMB2017* folder and execute the following commands to call footprints using ATAC-seq data:
::
    
    rgt-hint --atac-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=B ./session2/step1/input/B_ATAC_chr1.bam ./session2/step1/input/B_ATACPeaks_chr1.bed
    rgt-hint --atac-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=CD4 ./session2/step1/input/CD4_ATAC_chr1.bam ./session2/step1/input/CD4_ATACPeaks_chr1.bed
    rgt-hint --atac-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=LSK ./session2/step1/input/LSK_ATAC_chr1.bam ./session2/step1/input/LSK_ATACPeaks_chr1.bed

To call footprints with histone data execute the commands:
::
    rgt-hint --histone-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=B_H3K27Ac_chr1_footprints ./session2/step1/B_H3K27Ac_chr1.bam ./session2/step1/B_H3K27AcPeaks_chr1.bed
    rgt-hint --histone-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=CD4_H3K27Ac_chr1_footprints ./session2/step1/CD4_H3K27Ac_chr1.bam ./session2/step1/CD4_H3K27AcPeaks_chr1.bed
    rgt-hint --histone-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=LSK_H3K27Ac_chr1_footprints ./session2/step1/LSK_H3K27Ac_chr1.bam ./session2/step1/LSK_H3K27AcPeaks_chr1.bed


Step2: Intersecting footprints with differential histone peaks
--------------------------------------------------------------

To derive candidate regions for TF binding, we combine (1) genome wide footprint calls and (2) genome wide differential histone peak calls using
the active chromatin marks H3K4me3 and H3K27ac. 

Step3: Deriving candidate transcriptional regulators 
----------------------------------------------------


