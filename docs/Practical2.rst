==================================================================
Practical II - Footprint calling & Transcription factor prediction
==================================================================


Footprint calling
-----------------------------------------------

After installation of HINT, please go to EpigenomicsTutorial-ISMB2017 folder and execute the following commands to perform footprints identification using ATAC-seq
::
    
    rgt-hint --organism=mm10 --atac-footprints --output-location=./output/session2 --output-prefix=B ./input/session2/B_ATAC_chr1.bam ./input/session2/B_ATACPeaks_chr1.bed
    rgt-hint --organism=mm10 --atac-footprints --output-location=./output/session2 --output-prefix=CD4 ./input/session2/CD4_ATAC_chr1.bam ./input/session2/CD4_ATACPeaks_chr1.bed
    rgt-hint --organism=mm10 --atac-footprints --output-location=./output/session2 --output-prefix=LSK ./input/session2/LSK_ATAC_chr1.bam ./input/session2/LSK_ATACPeaks_chr1.bed

using histone data
::
    rgt-hint --organism=mm10 --histone-footprints --output-location=./output/session2 --output-prefix=B_H3K27Ac_chr1_footprints ./input/session2/B_H3K27Ac_chr1.bam ./input/session2/B_H3K27AcPeaks_chr1.bed
    rgt-hint --organism=mm10 --histone-footprints --output-location=./output/session2 --output-prefix=CD4_H3K27Ac_chr1_footprints ./input/session2/CD4_H3K27Ac_chr1.bam ./input/session2/CD4_H3K27AcPeaks_chr1.bed
    rgt-hint --organism=mm10 --histone-footprints --output-location=./output/session2 --output-prefix=LSK_H3K27Ac_chr1_footprints ./input/session2/LSK_H3K27Ac_chr1.bam ./input/session2/LSK_H3K27AcPeaks_chr1.bed

