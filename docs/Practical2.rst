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

To perform footprinting, we need (1) a bam file containing the sequence data and (2) a bed file including the interesting regions, such as DNase I hypersensitive sites (DHSs) for DNase-seq and peaks for ATAC-seq. 

All input files are found in the folder ``/EpigenomicsTutorial-ISMB2017/session2/step1/input``

**1.** Assure that you are in the directory ``EpigenomicsTutorial-ISMB2017/output/session2``, otherwise *cd* to that directory.

**2.** Generate an output folder for the resulting files and **enter the folder**:
::
    mkdir step1
    cd step1

**3.** Execute the following commands to call footprints using ATAC-seq data:
::
    rgt-hint --atac-footprints --organism=mm10 --output-location=./ --output-prefix=B_ATAC_chr1_footprints ../../../session2/step1/input/B_ATAC_chr1.bam ../../../session2/step1/input/B_ATACPeaks_chr1.bed
    rgt-hint --atac-footprints --organism=mm10 --output-location=./ --output-prefix=CD4_ATAC_chr1_footprints ../../../session2/step1/input/CD4_ATAC_chr1.bam ../../../session2/step1/input/CD4_ATACPeaks_chr1.bed
    rgt-hint --atac-footprints --organism=mm10 --output-location=./ --output-prefix=LSK_ATAC_chr1_footprints ../../../session2/step1/input/LSK_ATAC_chr1.bam ../../../session2/step1/input/LSK_ATACPeaks_chr1.bed

To call footprints with histone data execute the commands:
::
    rgt-hint --histone-footprints --organism=mm10 --output-location=./ --output-prefix=B_H3K27Ac_chr1_footprints ../../../session2/step1/input/B_H3K27Ac_chr1.bam ../../../session2/step1/input/B_H3K27AcPeaks_chr1.bed
    rgt-hint --histone-footprints --organism=mm10 --output-location=./ --output-prefix=CD4_H3K27Ac_chr1_footprints ../../../session2/step1/input/CD4_H3K27Ac_chr1.bam ../../../session2/step1/input/CD4_H3K27AcPeaks_chr1.bed
    rgt-hint --histone-footprints --organism=mm10 --output-location=./ --output-prefix=LSK_H3K27Ac_chr1_footprints ../../../session2/step1/input/LSK_H3K27Ac_chr1.bam ../../../session2/step1/input/LSK_H3K27AcPeaks_chr1.bed

The above commands will output a BED file containing the footprints, inside the current folder. Each footprint, i.e. each line of the BED file, will also contain information regarding the tag-count score of each footprint. This score can be used as a footprint quality assessment (the higher the value, the better). In addition, a file including the details of reads and footprints will also be written in the same folder of BED file.

**4.** Execute the following commands to do motif matching:
::
    rgt-motifanalysis --matching --organism=mm10 --output-location=./ --use-only-motifs=../../../session2/step1/input/motifs.txt ../../../session2/step1/result/B_ATAC_footprints.bed
    rgt-motifanalysis --matching --organism=mm10 --output-location=./ --use-only-motifs=../../../session2/step1/input/motifs.txt ../../../session2/step1/result/CD4_ATAC_footprints.bed
    rgt-motifanalysis --matching --organism=mm10 --output-location=./ --use-only-motifs=../../../session2/step1/input/motifs.txt ../../../session2/step1/result/Lsk_ATAC_footprints.bed

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

Please check the `documentation <https://github.com/SchulzLab/TEPIC/blob/master/docs/Description.pdf>`_ for details on the method.

We provide a script that automatically performs steps (1) and (2) as well as necessary data processing and formatting steps (See `DYNAMITE documentation <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_ for details).
All files used in this step are available in ``/EpigenomicsTutorial-ISMB2017/session2/Step3/input``. Additionally, we require the mm10 reference genome, which you should have downloaded while installing *HINT*.

Note that we precomputed the differential gene expression estimates. Computing those is neither part of the actual tutorial nor of the *DYNAMITE* workflow. However a tool you could use to compute differential gene/transcript expression is `Cuffdiff <http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/>`_.

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
The results of the analysis will be stored seperately for each run in ``EpigenomicsTutorial-ISMB2017/output/session2/step3/``. There are three subfolders for
each comparison:

#.Affinities
#.IntegratedData
#.Learning_Results

The folder *Affinities* contains TF affinities calculated in the provided regions for both groups, gene TF scores for both groups, and a metadata file that
lists all settings used for the TF annotation with *TEPIC* (subfolders *group1* and *group2*). The subfolder *mean* contains the mean gene TF scores computed
for group1 and group2. This is needed if you analyse more than one biological replicate per group. The folder *ratio* contains the gene TF score ratios computed between
the gene TF scores of group1 and group2.

The folder *IntegratedData* encloses matrices that are composed of (1) gene TF score ratios and (2) a measure of differential gene expression. In the folder *Log2* the differential gene expression
is represented as the log2 expression ratio between group1 and group2. In the folder *Binary* a 1 means a gene is upregulated in group 1 compared to group 2, whereas a 0 means it is downregulated in group1.
The binary format is used as input for the classification. 

The folder *Learning_Results* comprises the results of the logistic regression classifier. The following files should be produced if all R dependencies are available:

#. Performance_overview.txt
#. Confusion-Matrix_<1..6>_Integrated_Data_For_Classification.txt
#. Regression_Coefficients_Cross_Validation_Integrated_Data_For_Classification.txt
#. Regression_Coefficients_Entire_Data_Set_Integrated_Data_For_Classification.txt
#. Performance_Barplots.png
#. Regression_Coefficients_Cross_Validation_Heatmap_Integrated_Data_For_Classification.svg
#. Regression_Coefficients_Entire_Data_SetIntegrated_Data_For_Classification.png
#. Misclassification_Lambda_<1..6>_Integrated_Data_For_Classification.svg

The file *Performance_overview.txt* contains accuracy on Test and Training data sets as well as F1 measures. These values are visualised in *Performance_Barplots.png*.
As the name suggests, the files *Confusion-Matrix_<1..6>_Integrated_Data_For_Classification.txt* contain the confusion matrix computed on the test data sets.
They show model performance by reporting True Positives (TP), False Positives (FP), True Negatives (TN), and False Negatives (FN) in the following layout:

+---------------------+----------+----------+
| Observed/Predicted  | Positive | Negative |
+=====================+==========+==========+
| Positive            |    TP    |    FN    |
+---------------------+----------+----------+
| Negative            |    FP    |    TN    |
+---------------------+----------+----------+

The heatmap *Regression_Coefficients_Cross_Validation_Heatmap_Integrated_Data_For_Classification.svg* shows the regression coefficients of all selected features in
the outer cross validation. This is very well suited to find features that are stably selected in all outer cross validation folds. The raw data used to generate the figure is stored in 
*Regression_Coefficients_Cross_Validation_Integrated_Data_For_Classification.txt*. The stronger a regression coefficient, the more important it is in the model.

In addition to the heatmap showing the regression coefficients during the outer cross validation, we also compute the regression coefficients that are learned on the full
data set: *Regression_Coefficients_Entire_Data_SetIntegrated_Data_For_Classification.png* and *Regression_Coefficients_Entire_Data_Set_Integrated_Data_For_Classification.txt*.

The figures *Misclassification_Lambda_<1..6>_Integrated_Data_For_Classification.svg* are of technical nature. They show the relationship between the misclassification error and the lambda parameter of the logistic regression function. 

**5.** In addition to the plots describing model performance and feature selection generated by *DYNAMITE* (as described `here <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_), you can create further Figures for a distinct feature of interest
using the script ``TEPIC/MachineLearningPipelines/DYNAMITE/Scripts/generateFeaturePlots.R``. This will provide you with density plots showing the distribution of the feature in 
the two cell types, scatter plots linking feature value to gene expression changes, and a mini heatmap visualising the features regression coefficients. 

To use this script, go to the folder ``TEPIC/MachineLearningPipelines/DYNAMITE/Scripts/`` and use the command
::

	Rscript generateFeaturePlots.R /local/home/fschmidt/Documents/Research/EpigenomicsTutorial-ISMB2017/output/session2/step3/LSK-vs-CD4/ HOXA3 LSK CD4


This command will generate a plot comparing HOXA3 in LSK against CD4. Feel free to look at further features as you wish. The figure will be stored in the specified directory that contains the results of the *DYNAMITE* analysis.
Again, note that you have to **replace** the prefix ``/local/home/fschmidt/Documents/Research/`` with the proper path used on your system. 
Precomputed results are stored in ``/EpigenomicsTutorial-ISMB2017/session2/Step3/result`.
