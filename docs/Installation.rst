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
* `samtools <http://samtools.sourceforge.net>`_

**Note:** that the individual softwares may have some other dependencies, e.g. bedtools which you should have installed.

============
Installation
============

The following software packages need to be installed for running the tutorial:

R version 3.2 or higher.

`histoneHMM <https://github.com/matthiasheinig/histoneHMM>`_ 
-----------------------------------------------

You might need to install the following dependencies before installing histoneHMM.

:strong:`Unix libraries:`
  * lib-gcc

:strong:`R libraries:`
  * Rcpp
  * GenomicRanges (bioconductor.org)
  * Rsamtools (bioconductor.org)
  * mvtnorm

The package was developed and tested using a linux system and R. 
To install the latest version of the package, open an R terminal and type in the following commands (using the 'devtools' package):

**In the R terminal type the following**::

  install.packages("devtools") # if devtools is not yet installed
  devtools::install_github("matthiasheinig/histoneHMM")

Now the latest version of histoneHMM should be installed on your system.
In the tutorial, we will use the command-line interface to histoneHMM. In order for this to work smoothly, it would be best if you add the path to the histoneHMM script files to your $PATH variable (otherwise you'd have to specify the full path each time you call histoneHMM). You can find out the path you need to add by starting an R terminal and typing:
::
  system.file("bin/", package="histoneHMM")

Which should yield something like this:
::
  [1] "/home/[user-path]/R/x86_64-redhat-linux-gnu-library/3.2/histoneHMM/bin/"

Now you can add the directory indicated above to your PATH variable by calling:
::
  export PATH=$PATH:/home/[user-path]/R/x86_64-redhat-linux-gnu-library/3.2/histoneHMM/bin/

If you want to make the histoneHMM command-line available to you everytime you log on to your system, make sure that the directory is added to the $PATH variable everytime you log on or create a new terminal (e.g. by modifying your ~/.bashrc).


`HINT <http://github.com/CostaLab/reg-gen>`_ 
-----------------------------------------------

To install HINT (RGT Suite), you are advised to use the Python package installer pip. First, download the pip installer `get-pip.py <http://bootstrap.pypa.io/get-pip.py>`_ and then install pip.

::

    python get-pip.py

Next, install dependencies:

::

    pip install --user cython numpy scipy
    pip install --user https://github.com/fabioticconi/MOODS/tarball/pypi-ready

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
