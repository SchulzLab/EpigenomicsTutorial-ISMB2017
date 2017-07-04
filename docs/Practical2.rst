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
    rgt-hint --histone-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=B_H3K27Ac_chr1_footprints ./session2/step1/input/B_H3K27Ac_chr1.bam ./session2/step1/input/B_H3K27AcPeaks_chr1.bed
    rgt-hint --histone-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=CD4_H3K27Ac_chr1_footprints ./session2/step1/input/CD4_H3K27Ac_chr1.bam ./session2/step1/input/CD4_H3K27AcPeaks_chr1.bed
    rgt-hint --histone-footprints --organism=mm10 --output-location=./output/session2/step1 --output-prefix=LSK_H3K27Ac_chr1_footprints ./session2/step1/input/LSK_H3K27Ac_chr1.bam ./session2/step1/input/LSK_H3K27AcPeaks_chr1.bed


Step2: Intersecting footprints with differential histone peaks
-----------------------------------------------

To derive candidate regions for TF binding, we combine (1) genome wide footprint calls and (2) genome wide differential histone peak calls using
the active chromatin marks H3K4me3 and H3K27ac. In addition to default unix functions we  use *bedtools* to combine the respective bed files. 

All input files are available in the folder ``/EpigenomicsTutorial-ISMB2017/session2/Step2/input``.

1. Generate an output folder for the resulting bed files and enter the folder:
::

	cd EpigenomicsTutorial-ISMB2017/output/session2/
	mkdir step2
	cd step2

2. Combine the Differential peak calls for H3K4me3 and H3K27ac in one, sorted bed file. This needs to be done for each pairwise comparison and each cell type individually:
::

	cat ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K27ac-vs-CD4_H3K27ac-B.bed ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K4me3-vs-CD4_H3K4me3-B.bed | sort -k1,1 -k2,2n > B_vs_CD4_H3K27ac_H3K4me3_B_sorted.bed
	cat ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K27ac-vs-CD4_H3K27ac-CD4.bed ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K4me3-vs-CD4_H3K4me3-CD4.bed | sort -k1,1 -k2,2n > B_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed

	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-B_H3K27ac-LSK.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-B_H3K4me3-LSK.bed | sort -k1,1 -k2,2n > LSK_vs_B_H3K27ac_H3K4me3_LSK_sorted.bed
	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-B_H3K27ac-B.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-B_H3K4me3-B.bed | sort -k1,1 -k2,2n > LSK_vs_B_H3K27ac_H3K4me3_B_sorted.bed

	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-CD4_H3K27ac-LSK.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-CD4_H3K4me3-LSK.bed | sort -k1,1 -k2,2n > LSK_vs_CD4_H3K27ac_H3K4me3_LSK_sorted.bed
	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-CD4_H3K27ac-CD4.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-CD4_H3K4me3-CD4.bed | sort -k1,1 -k2,2n > LSK_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed


The *cat* command aggregates the input files for H3K27ac and H3K4me3 and pipes them (using the *|* operator) to a sort function which sorts by chromosome (k1,1) and first genomic corrdinate (k2,2n). 
The result is stored in a specified output bed file (using the *>* operator).

3. Merge overlapping histone peaks using *bedtools merge* and intersect the merged regions with HINT-BCs footprint calls using *bedtools intersect*:
::
	
	bedtools merge -i B_vs_CD4_H3K27ac_H3K4me3_B_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/B.bed > Footprints_B_vs_CD4_H3K27ac_H3K4me3_B.bed
	bedtools merge -i B_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/CD4.bed > Footprints_B_vs_CD4_H3K27ac_H3K4me3_CD4.bed

	bedtools merge -i LSK_vs_CD4_H3K27ac_H3K4me3_LSK_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/LSK.bed > Footprints_LSK_vs_CD4_H3K27ac_H3K4me3_LSK.bed
	bedtools merge -i LSK_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/CD4.bed > Footprints_LSK_vs_CD4_H3K27ac_H3K4me3_CD4.bed

	bedtools merge -i LSK_vs_B_H3K27ac_H3K4me3_LSK_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/LSK.bed > Footprints_LSK_vs_B_H3K27ac_H3K4me3_LSK.bed
	bedtools merge -i LSK_vs_B_H3K27ac_H3K4me3_B_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/B.bed > Footprints_LSK_vs_B_H3K27ac_H3K4me3_B.bed

The *bedtools merge* command combines to overlapping regions into one region. The result of the intersection is piped into the standard input stream (*stdin*) of the *bedtools intersect -a* argument, while the *-b* argument
is result of the Footprint calling. The resulting files will contain only footprints that intersect with a differential H3K4me3 and/or H3K27ac peak. In the next step, we will use these regions as candidate regions for TF binding. 

Step3: Deriving candidate transcriptional regulators 
----------------------------------------------------

