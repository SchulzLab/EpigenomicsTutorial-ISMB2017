==================================================================
Practical II - Footprint calling & Transcription factor prediction
==================================================================
In the second practical, we will perform a footprint analysis with `HINT <http://www.regulatory-genomics.org/hint/>`_ to identify cell specific binding sites from open chromatin data (ATAC-seq). Next, we 
will combine the footprints with the differential histone peaks detected by `histoneHMM <http://histonehmm.molgen.mpg.de>`_ (c.f. practical 1). 
Thereby, we will find tissue specific TF binding sites, which are located in regions with cell spefici histone peaks. These regulatory regions are used in a 
`DYNAMITE <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_ analysis with the aim
of inferring TFs might be related to gene expression differences between the tissues of interest. 

Step1: Footprint calling
-----------------------------------------------

First, we will use `HINT <http://www.regulatory-genomics.org/hint/>`_ to find genomic regions (footprints) with cell specific TF binding sites. For this, HINT requires (1) a sorted bam file containing the aligned reads from the sequencing library (DNase-,ATAC- or histone ChIP-seq) (2) and a bed file including peaks detected in the same sequencing library provided in (1). These peak regions are used by HINT to reduce the search space and can be geneated by any  peak caller. 

Here, we will analyse ATAC-seq data from LSK cells (equivalent to MPP cells), B cells and T CD4 cells obtained from `Lara-Astiaso et al 2014 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60103>`_. We have prepared low level analysis files including read alignment and peaks calling on chromossome 1, which can be found in this folder ``/EpigenomicsTutorial-ISMB2017/session2/step1/input``. Check here for an script describing how these were generated ``XXX - to do''.

**1.** First, go to the EpigenomicsTutorial-ISMB2017 directory and generate an output folder for the result files:
::
    cd EpigenomicsTutorial-ISMB2017
    mkdir session2/step1/output
   
**2.** Execute the following commands to call footprints on B, LSK and CD4 cells:
::
    rgt-hint --atac-footprints --organism=mm10 --output-location=session2/step1/output/ --output-prefix=B_ATAC_chr1_footprints session2/step1/input/B_ATAC_chr1.bam session2/step1/input/B_ATACPeaks_chr1.bed
    rgt-hint --atac-footprints --organism=mm10 --output-location=session2/step1/output/ --output-prefix=LSK_ATAC_chr1_footprints session2/step1/input/LSK_ATAC_chr1.bam session2/step1/input/LSK_ATACPeaks_chr1.bed
    rgt-hint --atac-footprints --organism=mm10 --output-location=session2/step1/output/ --output-prefix=CD4_ATAC_chr1_footprints session2/step1/input/CD4_ATAC_chr1.bam session2/step1/input/CD4_ATACPeaks_chr1.bed

This will generate an output file  ``session2/step1/output/B_ATAC_chr1_footprints.bed`` containing the genomic locations of the footprints. HINT also produces a file with ending ".info", which has general statistics from the analysis as no. of footprints, total number of reads and so on. 

We can use IGV to vizualise ATAC-seq signals and footprint predictions in particular loci. First, we can use a special HINT command to generate genomic profiles (bigWig files). 
::
   rgt-hint XXX 

This bigwig file contains number of ATAC-seq (or DNAse-seq) reads at each genomic position as estimated by HINT after signal normalization and cleveage bias correction. This is therefore more accurate than simply looking a coverage profiles of a bam file. Open all bw and bam files from ATAC-seq in your IGV. Remember to set the genome version to mm10 beforehand. You can also enrich the data by opening bam files of histone modifications as H327ac of the same cells. Go to the gene Ilr1 and zoom in on its promoter. We observe that this gene only has signals in LSK cells. Moreover, ATAC-seq reads only covers a small part of the H327ac, which is acessible for binding. (XXX - choose relevant gene and include a screenshot here). 

**3.** A import question when doing footprint analysis is to evaluate which Tf motifs overllap with footprints and evaluate the ATAC-seq profiles around these motifs. RGT suite also offers a tool for finding motif matches. For example, we analyse here motifs from factors SPI1 and ELK4, which were discussed in Lara-Astiaso et al. 2014 to be associated respectivelly associated to LSK and B cells.

Execute the following commands to do motif matching inside footprints:
::
    rgt-motifanalysis --matching --organism=mm10 --output-location=session2/step1/output/ --use-only-motifs=session2/step1/input/motifs.txt session2/step1/result/B_ATAC_footprints.bed
    rgt-motifanalysis --matching --organism=mm10 --output-location=session2/step1/output/ --use-only-motifs=session2/step1/input/motifs.txt session2/step1/result/CD4_ATAC_footprints.bed
    rgt-motifanalysis --matching --organism=mm10 --output-location=session2/step1/output/ --use-only-motifs=session2/step1/input/motifs.txt session2/step1/result/Lsk_ATAC_footprints.bed

The file ``session2/step1/motifs.txt``  contains a list of `JASPAR <http://jaspar.genereg.net/>`_ motif ids to be used in the analysis. Ignoring this option will search for all JASPAR motifs. The above commands will generate three BED files (i.e. LSK_ATAC_footprints_mpbs.bed) containing the matched motif instances overllaping with distinct footprint regions. The 4th collumn contains the motif name and the 5th collumn the bitscore (higher values indicates best match).  You can open MPBS.bed files in IGV and you will observe that there is one SFP1 binding site overllaping on the promoter of Ilr1.

**4.** Finally, we use HINT to generate average ATAC-seq profiles around binding sites of particular TF. These analysis allow us to inspect the cut profiles and the underlying sequence conservation. Moreover, by comparing the cut profiles from two ATAC-seq libraries (i.s. LSK vs B cells), it is possible to inspect the chromatin and TF activity status. For this, execute the following commmands:
::
    mkdir session2/step1/output/LSK_B
    rgt-hint --diff-footprints --organism=mm10 --mpbs-file=session2/step1/output/B_ATAC_footprints_mpbs.bed --reads-file1=session2/step1/input/B_ATAC.bam --reads-file2=session2/step1/input/LSK_ATAC.bam --output-location=session2/step1/output/LSK_B --output-prefix=LSK_B

    mkdir session2/step1/output/LSK_CD4
    rgt-hint --diff-footprints --organism=mm10 --mpbs-file=session2/step1/output/CD4_ATAC_footprints_mpbs.bed --reads-file1=session2/step1/input/B_ATAC.bam --reads-file2=session2/step1/input/LSK_ATAC.bam --output-location=session2/step1/output/LSK_CD4 --output-prefix=LSK_CD4

The above commands will generate pdf (and eps) files with a ATAC-seq profile for each of the motifs founds in the input bed files. Let's check the profiles in the comparirson LSK and CD4, you will see that ELK4 has higher number of ATAC-seq counts in CD4 cells, while SFP1 has more ATAC-seq in LSK cells. This fits with the results discussed in Lara-Astiaso that SFP1 are more relevant/active in LSK, while ELK4 in CD4 cells.

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
Precomputed results are stored in ``/EpigenomicsTutorial-ISMB2017/session2/Step3/result``.
