=============
What you need
=============

You need to bring your own **laptop** with the following software installed (see detailed instructions below)

* R version 3.2 or higher
* Python version 2.7
* `histoneHMM <http://histonehmm.molgen.mpg.de>`_ 
* `HINT <http://github.com/CostaLab/reg-gen>`_ 
* `TEPIC <https://github.com/SchulzLab/TEPIC>`_ 
* `IGV <http://software.broadinstitute.org/software/igv/>`_

**Note:** that the individual softwares may have some other dependencies, e.g. bedtools which you should have installed.

============
Installation
============

The following software packages need to be installed for running the tutorial:

R version 3.2 or higher.

`histoneHMM <http://histonehmm.molgen.mpg.de>`_ 
-----------------------------------------------

You might need to install the following dependencies before installing histoneHMM.

:strong:`Unix libraries:`
  * lib-gcc

:strong:`R libraries:`
  * Rcpp
  * GenomicRanges (bioconductor.org)
  * Rsamtools (bioconductor.org)
  * mvtnorm

The package was developed and tested using a linux system, so installation instructions are given for linux and unix like systems.

**Start a terminal and download as follows**::

  wget http://histonehmm.molgen.mpg.de/histoneHMM_1.6.tar.gz


**Install using**::

  R CMD INSTALL histoneHMM_1.6.tar.gz

`HINT <http://github.com/CostaLab/reg-gen>`_ 
-----------------------------------------------

To install HINT (RGT Suite), you are advised to use the Python package installer pip. First, download the pip installer `get-pip.py <http://bootstrap.pypa.io/get-pip.py>`_ and then install pip.

::

    python get-pip.py

Next, install dependencies:

::

    pip install --user cython numpy scipy


and finally install HINT and RGT suite.

::

    pip install --user RGT

Alternatively, look at detailed instructions `here <http://www.regulatory-genomics.org/hint/introduction/>`_.

You also need to download genome information for mouse genome mm10.

::

    cd ~/rgtdata
    python setupGenomicData.py --mm10


`IGV <http://software.broadinstitute.org/software/igv/>`_
-----------------------------------------------

Instructions on installing IGV are available `here <http://software.broadinstitute.org/software/igv/download>`_. We advise you to download a binary distribution. 

`TEPIC <https://github.com/SchulzLab/TEPIC>`_ 
-----------------------------------------------

**Dependencies**

TEPIC requires:

  * bedtools
  * A C++ compiler supporting openmp, e.g. g++ (test with version 4.9.2)
  
To run the machine learning pipeline DYNAMITE, which is part of the TEPIC repository, we require the `R libraries:`

  * glmnet
  * doMC
  * gplots
  * ggplot2
  * reshape2
  * gridExta
  
The TEPIC examples in the tutorial also require the mouse reference genome that was downloaded during the HINT setup. 

**Installation**

Start a terminal and clone the TEPIC repository ::

  git clone https://github.com/SchulzLab/TEPIC.git
  
Next, go to the folder ::

  TEPIC/Code
  
and type ::

  bash compileTRAP.sh
  
to build the C++ component of TEPIC.

If all dependencies mentioned above are available, no further installation steps are required. 

**Testing**

To test the core functionality of TEPIC, go to the folder::
   
   TEPIC/Code/ 
   
and run the example with the command:::

  ./TEPIC.sh -g ../Test/example_sequence.fa -b ../Test/example_regions.bed -o TEPIC-Example -p ../PWMs/pwm_vertebrates_jaspar_uniprobe_original.PSEM -a ../Test/example_annotation.gtf -w 3000 -e FALSE

There should be three result files generated:

  * TEPIC-Example <date> Affinity.txt
  * TEPIC-Example <date> amd.tsv
  * TEPIC-Example <date> Peak_Features_Affinity_Gene_View_Filtered.txt
  
To test the logistic regression framework DYNAMITE, which will be used in the tutorial, go to the folder ::

  /TEPIC/MachineLearningPipelines/DYNAMITE/
  
and run the provided example by entering the command ::

  bash runDYNAMITE.sh ./DYNAMITE.cfg
  
This will generate all output files that are described in the `DYNAMITE documentation <https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/README.md>`_. 

For further information, please see the `TEPIC repository <https://github.com/SchulzLab/TEPIC>`_ . 
