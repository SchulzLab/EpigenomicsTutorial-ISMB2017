==================================================================
Practical II - Footprint calling & Transcription factor prediction
==================================================================
In the second practical, we will first identify TF-footprints using the software `HINT <http://github.com/CostaLab/reg-gen>`_ and
combine those with differential histone peaks called using `histoneHMM <http://histonehmm.molgen.mpg.de>`_ (c.f. practical 1).
Thereby, we will identify tissue specific sets of differential candidate binding regions for TF binding. These are used in a 
`DYNAMITE documentation <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_ analysis with the aim
of inferring TFs might be related to gene expression differences between the tissues of interest. 

Step1: Footprint calling
-----------------------------------------------

After installation of HINT, go to the ``EpigenomicsTutorial-ISMB2017`` folder and execute the following commands to call footprints using ATAC-seq data:
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

All input files are available in the folder ``/EpigenomicsTutorial-ISMB2017/session2/step2/input``.

**1.** Assure that you are in the directory ``EpigenomicsTutorial-ISMB2017/output/session2``, otherwise *cd* to that directory.

**2.** Generate an output folder for the resulting bed files and **enter the folder**:
::
	mkdir step2
	cd step2
	
**3.** Combine the Differential peak calls for H3K4me3 and H3K27ac in one, sorted bed file. This needs to be done for each pairwise comparison and each cell type individually:
::
	cat ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K27ac-vs-CD4_H3K27ac-B.bed ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K4me3-vs-CD4_H3K4me3-B.bed | sort -k1,1 -k2,2n > B_vs_CD4_H3K27ac_H3K4me3_B_sorted.bed
	cat ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K27ac-vs-CD4_H3K27ac-CD4.bed ../../../session2/step2/input/Dif_Histone_Peaks/B_H3K4me3-vs-CD4_H3K4me3-CD4.bed | sort -k1,1 -k2,2n > B_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed

	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-B_H3K27ac-LSK.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-B_H3K4me3-LSK.bed | sort -k1,1 -k2,2n > LSK_vs_B_H3K27ac_H3K4me3_LSK_sorted.bed
	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-B_H3K27ac-B.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-B_H3K4me3-B.bed | sort -k1,1 -k2,2n > LSK_vs_B_H3K27ac_H3K4me3_B_sorted.bed

	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-CD4_H3K27ac-LSK.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-CD4_H3K4me3-LSK.bed | sort -k1,1 -k2,2n > LSK_vs_CD4_H3K27ac_H3K4me3_LSK_sorted.bed
	cat ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K27ac-vs-CD4_H3K27ac-CD4.bed ../../../session2/step2/input/Dif_Histone_Peaks/LSK_H3K4me3-vs-CD4_H3K4me3-CD4.bed | sort -k1,1 -k2,2n > LSK_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed

The *cat* command aggregates the input files for H3K27ac and H3K4me3 and pipes them (using the *|* operator) to a sort function which sorts by chromosome (*k1,1*) and first genomic coordinate (*k2,2n*). The result is stored in a specified output bed file (using the *>* operator).

**4.** Merge overlapping histone peaks using *bedtools merge* and intersect the merged regions with HINT-BCs footprint calls using *bedtools intersect*:
::
	
	bedtools merge -i B_vs_CD4_H3K27ac_H3K4me3_B_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/B.bed > Footprints_B_vs_CD4_H3K27ac_H3K4me3_B.bed
	bedtools merge -i B_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/CD4.bed > Footprints_B_vs_CD4_H3K27ac_H3K4me3_CD4.bed

	bedtools merge -i LSK_vs_CD4_H3K27ac_H3K4me3_LSK_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/LSK.bed > Footprints_LSK_vs_CD4_H3K27ac_H3K4me3_LSK.bed
	bedtools merge -i LSK_vs_CD4_H3K27ac_H3K4me3_CD4_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/CD4.bed > Footprints_LSK_vs_CD4_H3K27ac_H3K4me3_CD4.bed

	bedtools merge -i LSK_vs_B_H3K27ac_H3K4me3_LSK_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/LSK.bed > Footprints_LSK_vs_B_H3K27ac_H3K4me3_LSK.bed
	bedtools merge -i LSK_vs_B_H3K27ac_H3K4me3_B_sorted.bed | bedtools intersect -a stdin -b ../../../session2/step2/input/Footprints/B.bed > Footprints_LSK_vs_B_H3K27ac_H3K4me3_B.bed

The *bedtools merge* command combines to overlapping regions into one region. The result of the intersection is piped into the standard input stream (*stdin*) of the *bedtools intersect -a* argument, while the *-b* argument
is result of the Footprint calling. The resulting files will contain only footprints that intersect with a differential H3K4me3 and/or H3K27ac peak. In the next step, we will use these regions as candidate regions for TF binding. 
Precomputed results are stored in ``/EpigenomicsTutorial-ISMB2017/session2/Step2/result``.


Step3: Deriving candidate transcriptional regulators using *DYNAMITE*
----------------------------------------------------

During a *DYNAMITE* analysis, two main computational tasks are undertaken:

#. We calculate TF binding affinities for an example data set of 93 TFs and aggregate those to gene-TF scores using *TEPIC*. TF affinities are a quantitative measure of TF binding to a distinct genomic region. 
#. A logistic regression classifier is learned that uses changes in TF gene scores between two samples to predict which genes are up/down- regulated between them. Investigating the features of the model allows the inference of potentially interesting regulators that are correlated to the observed expression changes. 

We provide a script that automatically performs steps (1) and (2) as well as necessary data processing and formatting steps (See `DYNAMITE documentation <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_ for details).
All files used in this step are available in ``/EpigenomicsTutorial-ISMB2017/session2/Step3/input``. Additionally, we require the mm10 reference genome, which you should have downloaded while installing *HINT*.

**1.** Assure that you are in the directory ``EpigenomicsTutorial-ISMB2017/output/session2``, otherwise *cd* to that directory.

**2.** Generate an output folder for the resulting files:
::
	mkdir step3
	
**3.** To run the *DYNAMITE* script go to the *DYNAMITE* folder in the *TEPIC* repository ``TEPIC/MachineLearningPipelines/DYNAMITE``.

**4.** Run the individual pairwise comparisons for LSK vs B
::
	
	bash runDYNAMITE.sh /local/home/fschmidt/Documents/Research/EpigenomicsTutorial-ISMB2017/session2/step3/input/DYNAMITE-LSKvsB.cfg

LSK vs CD4
::
	bash runDYNAMITE.sh /local/home/fschmidt/Documents/Research/EpigenomicsTutorial-ISMB2017/session2/step3/input/DYNAMITE-LSKvsCD4.cfg

and B vs CD4
::
	bash runDYNAMITE.sh /local/home/fschmidt/Documents/Research/EpigenomicsTutorial-ISMB2017/session2/step3/input/DYNAMITE-BvsCD4.cfg

Note that you have to **replace** the prefix ``/local/home/fschmidt/Documents/Research/`` with the proper path used on your system. 
The *cfg* files are configuration files that specify the path to all files needed in a *DYNAMITE* analysis, e.g. bed files for candidate binding regions.
The results of the analysis will be stored seperately for each run in ``EpigenomicsTutorial-ISMB2017/output/session2/step3/``.

**5.** In addition to the plots describing model performance and feature selection generated by *DYNAMITE* (as described `here <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_), you can create further Figures for a distinct feature of interest
using the script ``TEPIC/MachineLearningPipelines/DYNAMITE/Scripts/generateFeaturePlots.R``. This will provide you with density plots showing the distribution of the feature in 
the two cell types, scatter plots linking feature value to gene expression changes, and a mini heatmap visualising the features regression coefficients. 

To use this script, go to the folder ``TEPIC/MachineLearningPipelines/DYNAMITE/Scripts/`` and use the command
::

	Rscript generateFeaturePlots.R /local/home/fschmidt/Documents/Research/EpigenomicsTutorial-ISMB2017/output/session2/step3/LSK-vs-CD4/ HOXA3 LSK CD4


This command will generate a plot comparing HOXA3 in LSK vs CD4. Feel free to look at further features as you wish. The figure will be stored in the specified directory that contains the results of the *DYNAMITE* analysis.
Again, note that you have to **replace** the prefix ``/local/home/fschmidt/Documents/Research/`` with the proper path used on your system. 
Precomputed results are stored in ``/EpigenomicsTutorial-ISMB2017/session2/Step3/result`.
