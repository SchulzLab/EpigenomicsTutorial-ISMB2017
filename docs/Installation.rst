=============
What you need
=============

You need to bring your own **laptop** with the following software installed (see detailed instructions below)

* R version xxx or higher
* Python version 2.7
* `histoneHMM <http://histonehmm.molgen.mpg.de>`_ 
* `HINT </http://github.com/CostaLab/reg-gen>`_ 
* TEPIC
* IGV

**Note:** that the individual softwares may have some other dependencies, e.g. bedtools which you should have installed.

============
Installation
============

The following software packages need to be installed for running the tutorial:

R version xxx or higher.

`histoneHMM <http://histonehmm.molgen.mpg.de>`_ 

The package was developed and tested using a linux system, so installation instructions are given for linux and unix like systems.

**Start a terminal and download as follows**::

  wget http://histonehmm.molgen.mpg.de/histoneHMM_1.6.tar.gz


**Install using**::

  R CMD INSTALL histoneHMM_1.6.tar.gz

To install HINT/RGT, you are advised to use the Python package installer pip. First, download the pip installer `get-pip.py <https://bootstrap.pypa.io/get-pip.py>`

::

    python get-pip.py

Next, install dependencies:

::

    pip install --user cython numpy scipy


and then HINT/RGT. 

::

    pip install RGT

This will install the HINT and RGT suite with all dependencies. Alternatively, look at detailed instructions `here <http://www.regulatory-genomics.org/hint/download-installation/>`.




