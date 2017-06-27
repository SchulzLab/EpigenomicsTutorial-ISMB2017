=============
What you need
=============

You need to bring your own **laptop** with the following software installed (see detailed instructions below)

* R version xxx or higher
* Python version 2.7
* `histoneHMM <http://histonehmm.molgen.mpg.de>`_ 
* `HINT <http://github.com/CostaLab/reg-gen>`_ 
* TEPIC
* `IGV <http://software.broadinstitute.org/software/igv/>`_

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

`HINT <http://github.com/CostaLab/reg-gen>`_ 

To install HINT/RGT, you are advised to use the Python package installer pip. First, download the pip installer `get-pip.py <http://bootstrap.pypa.io/get-pip.py>`_ and then install pip.

::

    python get-pip.py

Next, install dependencies:

::

    pip install --user cython numpy scipy


and finally install HINT and RGT suite.

::

    pip install --user RGT

Alternatively, look at detailed instructions `here <http://www.regulatory-genomics.org/hint/download-installation/>`_.


`IGV <http://software.broadinstitute.org/software/igv/>`_

IGV installation instrucitons are found `here <http://software.broadinstitute.org/software/igv/download>`_. We advise you to download binary distribution. 



