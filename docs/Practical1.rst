===============================================
Practical I - Differential Histone peak calling
===============================================

In the first part of the practical, we will have a look at histone modifications in different cell-lines as measured by ChIP-seq experiments. We will use `histoneHMM <https://github.com/matthiasheinig/histoneHMM>`_ for calling regions in the genome which show histone modifications as well as for identifying those regions, which show differential modification states between cell-lines.
The data used in this part of the practical can be found in your checked out tutorial directory under ``/EpigenomicsTutorial-ISMB2017/session1/step1/input``

Step 1: Checking read alignments
-----------------------------------------------
Before we look at any modifications patterns in our ChIP-seq experiments, we shall get an impression on how our sequencing data looks like. 

**1.** First, change into the EpigenomicsTutorial-ISMB2017 directory and see which files are available as an input
::
  cd EpigenomicsTutorial-ISMB2017
  ls -lh session1/step1/input
  
You can see the *.bam and *.bai files for the three cell-lines and the two examined histone modifications. Now we want to get a brief overview on the nature of the bam files.

**2.** Create summary statistics using samtools
::
  mkdir session1/step1/output
  for i in session1/step1/input/*.bam ; do 
    samtools flagstat $i > session1/step1/output/$(basename $i .bam).summary ; 
  done

This will create for each of the available *.bam files a short read summary in the step1/output directory.
In those summary files you can see that the input files are rather low input. In our case, howvever, we can never the less work with those experiments rather well. 

**3.** Have a look at the *.bam files using the IGV

Just open the IGV, then via ``File->Load from File`` open your *.bam file of choice. Make sure that the Mouse (mm10) genome is selected in the upper left view of the browser, since this is the genome build which was used during read mapping. Examine the bam-tracks, what do you observe? Are there regions of high/low coverage?

Step 2: Calling modified regions
-----------------------------------------------
Now we know what we are dealing with and we are ready to begin the process of analyzing our histone marks using histoneHMM. The first step is to identify those regions in the genome which show histone modifications, i.e. to 'call regions'. histoneHMM works by first binning the genome into equally sized, non overlapping bins and then analyzing the number of reads falling into each of those bins. Necessary inputs for this step are 1) a chromosome lengths file, indicating length and name of chromosomes which are available and 2) the *.bam file from the ChIP-seq experiment of interest. Chromosome lengths files can easily be downloaded from 'UCSC goldenpath'. Alternatively, you can extract the information from the *.bam files directly, since it should always be encoded in the header (e.g. using ``samtools view -H test.bam``). Here, we will download the mm10 file from UCSCS, unzip it, and then run the command-line version of histoneHMM to call the modified regions. We also use the tool's -b parameter to set the size of the bins in which the genome should be devided to 2000bp.

NOTE: Before going on, make sure that the histoneHMM 'bin' directory is contained in you PATH variable (see installation instructions)

**1.** Run histoneHMM's 'call_regions'
::
  mkdir session1/step1/output/regions
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz | gunzip -c
  for i in session1/step1/input/*.bam ; do 
    prefix=session1/step1/output/regions/$(basename $i .bam)
    histoneHMM_call_regions.R -b 2000 -c chromInfo.txt -o ${prefix} $i &> ${prefix}.debug
  done

Now for each experiment, the script generated a set of files. Figure out what the different files are using the histoneHMM `manual <http://histonehmm.molgen.mpg.de/v1.6/histoneHMM.pdf>`_ . 
histoneHMM fits a mixture model to the counts using an EM algorithm. The two components of the mixture reflect two parts of the histogram: one with very high signal (high counts) and one with low signal values (low counts). Now check the generated count histograms, do you observe the two parts of the mixture fit? How does the count histogram look, would you have expected something like this?

Step 3: Differential region calling
-----------------------------------------------
more to come soon
