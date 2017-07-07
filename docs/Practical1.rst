===============================================
Practical I - Differential Histone peak calling
===============================================

In the first part of the practical, we will have a look at histone modifications in different cell-types as measured by ChIP-seq experiments in B-cell, CD4-cell and LSK (MPP) cell data from `Lara-Astiaso et al 2014 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60103>`_. Specifically, we look at the H3K4me3 and H3K27ac histone modifications and will analyze how they change between blood progenitor and more differentiated cells.
We will use `histoneHMM <https://github.com/matthiasheinig/histoneHMM>`_ for calling regions in the genome which show histone modifications as well as for identifying those regions, which show differential modification states between cell-types.
The data used in this part of the practical can be found in your checked out tutorial directory under ``/EpigenomicsTutorial-ISMB2017/session1/step1/input``

**The final version of the practical will be available at 19.07.2017 at the latest.**

Step 1: Checking read alignments
-----------------------------------------------
Before we look at any modifications patterns in our ChIP-seq experiments, we shall get an impression on how our sequencing data looks like. 

NOTE: In the step1 input directory, we also provide experiment files for the H3k4me1 and H3K4me2 histone marks. Those will not be fully processed using the scripts on this page, but you can look at them if you have any spare time left.

**1.** First, change into the EpigenomicsTutorial-ISMB2017 directory and see which files are available as an input
::
  cd EpigenomicsTutorial-ISMB2017/session1
  ls -lh step1/input
  
You can see the *.bam and *.bai files for the three cell-lines and the two examined histone modifications. You also see some *.wig files which we'll use later when looking at our data using IGV, ignore them for now. 
Now we want to get a brief overview on the nature of the bam files.

**2.** Create summary statistics using samtools
::
  mkdir -p step1/output/stats
  for i in step1/input/*.bam ; do 
    samtools flagstat $i > step1/output/stats/$(basename $i .bam).summary ; 
  done

This will create for each of the available *.bam files a short read summary in the step1/output/stats directory. 
Now check those files, what do you see? Also have a look at the header of the *.bam files, what can you observe?

**3.** Have a look at the *.bam files using the IGV

Just open the IGV, then via ``File->Load from File`` open your *.bam file of choice and the corresponding *.wig file. Make sure that the correct Mouse genome (mm10) is selected in the upper left view of the browser, since this is the genome build which was used during read mapping. Examine the loaded tracks, what do you observe? Are there regions of high/low coverage?

Step 2: Calling modified regions
-----------------------------------------------
Now we know what we are dealing with and we are ready to begin the process of analyzing our histone marks using histoneHMM. The first step is to identify those regions in the genome which show histone modifications, i.e. to 'call regions'. histoneHMM works by first binning the genome into equally sized, non overlapping bins and then analyzing the number of reads falling into each of those bins. Necessary inputs for this step are 1) a chromosome lengths file, indicating length and name of chromosomes which are available and 2) the *.bam file from the ChIP-seq experiment of interest. Chromosome lengths files can easily be downloaded from 'UCSC goldenpath'. Alternatively, you can extract the information from the *.bam files directly, since it should always be encoded in the header. You can use either of the two following code sections to get your chromosome lengths file.

**1.1.** Getting chromosome lengths from UCSC
::
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz
  gunzip chromInfo.txt.gz
  # we filter the chr1 only, since we only have chr1 reads
  grep chr1 chromInfo.txt > chromInfo.chr1.txt

**1.2.** Extracting chromosome lengths from *.bam files
::
  samtools view -H step2/input/B_H3k27ac.bam | grep SN:chr1 | cut -f 2,3 | sed s/[SL][NQ]://g > chromInfo.chr1.txt
  
With the chromosome lengths file in place, we now run the command-line version of histoneHMM to call the modified regions. We also use the tool's -b parameter to set the size of the bins in which the genome should be devided to 2000bp.

NOTE: Before going on, make sure that the histoneHMM 'bin' directory is contained in you PATH variable (see installation instructions)

**2.** Run histoneHMM's 'call_regions'
::
  mkdir -p step2/output/regions
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz
  gunzip chromInfo.txt.gz
  # we filter the chr1 only, since we only have chr1 reads
  grep chr1 chromInfo.txt > chromInfo.chr1.txt
  for i in step2/input/*.bam ; do 
    prefix=step2/output/regions/$(basename $i .bam)
    histoneHMM_call_regions.R -b 2000 -c chromInfo.chr1.txt -o ${prefix} $i &> ${prefix}.debug
  done

Now for each experiment, the script generated a set of files. Figure out what the different files are using the histoneHMM `manual <http://histonehmm.molgen.mpg.de/v1.6/histoneHMM.pdf>`_ . 
histoneHMM fits a mixture model to the counts using an EM algorithm. The two components of the mixture reflect two parts of the histogram: one with very high signal (high counts) and one with low signal values (low counts). Now check the generated count histograms, do you observe the two parts of the mixture fit? How does the count histogram look, would you have expected something like this?

Step 3: Differential region calling
-----------------------------------------------
The next and last step in this pipeline is formed by the differential region calling. Here we will compare experiments of the same histone modification in different cell-lines. 
To perform the differential region calling with histoneHMM, we only need a file with binned count information as is created during the previous step for both experiments we want to compare. 

NOTE: If you want you can redirect all output of histoneHMM using the '&>' operator as we did in the previous step.

**1.** Call differential regions
::
  odir=step3/output/differential
  mkdir -p ${odir}
  idir=step3/input/regions/
  
  # call differential analysis for all possible comparisons
  # for H3K4me3
  histoneHMM_call_differential.R --sample1 LSK_H3K4me3 --sample2 CD4_H3K4me3 --outdir ${odir} ${idir}/LSK_H3K4me3.txt ${idir}/CD4_H3K4me3.txt
  histoneHMM_call_differential.R --sample1 CD4_H3K4me3 --sample2 B_H3K4me3 --outdir ${odir} ${idir}/CD4_H3K4me3.txt ${idir}/B_H3K4me3.txt
  histoneHMM_call_differential.R --sample1 LSK_H3K4me3 --sample2 B_H3K4me3 --outdir ${odir} ${idir}/LSK_H3K4me3.txt ${idir}/B_H3K4me3.txt
  
  # for H3K27ac
  histoneHMM_call_differential.R --sample1 LSK_H3K27ac --sample2 CD4_H3K27ac --outdir ${odir} ${idir}/LSK_H3K27ac.txt ${idir}/CD4_H3K27ac.txt
  histoneHMM_call_differential.R --sample1 CD4_H3K27ac --sample2 B_H3K27ac --outdir ${odir} ${idir}/CD4_H3K27ac.txt ${idir}/B_H3K27ac.txt
  histoneHMM_call_differential.R --sample1 LSK_H3K27ac --sample2 B_H3K27ac --outdir ${odir} ${idir}/LSK_H3K27ac.txt ${idir}/B_H3K27ac.txt
  
histoneHMM again creates several output files (check the `manual <http://histonehmm.molgen.mpg.de/v1.6/histoneHMM.pdf>`_ do get to know those files). The infividual *.gff files contain the regions which are modified in both, none or only one of the compared experiments. For further analysis, we will only consider those regions which show an average posterior probability of at least 0.8. Also we want to make the *.gff files somewhat more convenient to deal with and convert them into *.bed files. You can do this however you want, here we will use a straight forward method using only Unix commands.

**2.** Filter and convert differential calls
::
  for i in step3/output/differential/*.gff ; do
    ofile=$(dirname $i)/$(basename $i .gff).post_08.bed
    awk '{split($9,arr,";"); split(arr[1],arr2,"="); }{if(arr2[2]>=0.8) print $1 "\t" $4-1 "\t" $5}' ${i} > ${ofile}
  done

The new *.bed files (with the .post_08 suffix) now contain the coordinates of the differential and modified/not modified regions for the analyzed experiment. To further get to know the results, check how many differential regions were discovered for each comparison after filtering. How many regions do you observe? Do the numbers differ between the individual histone marks?
As a last step, open again IGV and load the *.bam files as before. But now also add a few of the filtered *.bed files to add tracks which show e.g. the location of the differential peaks. Can you visually discern the differential peaks in the *.bam tracks? Do you agree with the results from histoneHMM?




