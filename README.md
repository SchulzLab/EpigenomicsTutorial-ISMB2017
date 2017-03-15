# Prediction of Regulatory Networks from Expression and Chromatin Data
## ISMB 2017 Tutorial

This is the accompanying website of the ISMB tutorial about epigenomics data analysis held in *Prague on the 21st of July*. It contains the schedule as well as information, and will be extended to contain all the data and software necessary for the tutorial.

# Organizers / Presenters

* [Dr. Ivan G. Costa](http://costalab.org/team-2/ivan-g-costa-group-leader-2/), RWTH Aachen University, Germany
* [Dr. Marcel Schulz](https://bioinf.mpi-inf.mpg.de/homepage/index.php?&account=mschulz), Saarland University & Max Planck Institute for Informatics, Germany
* [Dr. Matthias Heinig](https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/4158/Dr.-Heinig/index.html), Helmholtz Center Munich, Germany 


# Motivation
One of the main molecular mechanisms controlling the temporal and spatial 
expression of genes is transcriptional regulation. In this process, 
transcription factors (TFs) bind to the promoter and enhancers in the 
vicinity of a gene to recruit (or block) the transcriptional machinery 
and start gene expression. Inference of gene regulatory networks, i.e. 
factors controlling the expression of a particular gene, is a key challenge 
when studying development and disease progression. The availability of 
different experimental assays (Histone ChIP-seq, Dnase1-seq, ATAC-seq, 
NOME-seq etc.) that allow to map in-vivo chromatin dynamics and gene 
expression (RNA-seq), has triggered the development of novel computational 
modelling approaches for accurate prediction of TF binding and activity by 
integrating these diverse epigenomic datasets. However, in practice researchers 
are faced with the problems that come with handling diverse assays, 
understanding the tools involved and building specific workflows 
that are tailored to the data they have.

# Goals & Audience

This tutorial is targeted to an audience of bioinformaticians with previous 
experience in gene expression and next generation sequencing analysis. 
This Intermediary level tutorial will provide you knowledge on the use of 
state-of-art tools for inference of gene regulatory networks from chromatin 
and expression data. First, we will We will first review toolsthe state-of-the-art 
to 1) predict regulatory regions from different epigenetic datasets, e.g., 
using differential peak callers (histoneHMM - Heinig et al., 2015) or 
footprint methods (HINT - Gusmao et al., 2014) and 2) show how to determine 
cell-specific TF binding in these regions (PIQ - Sherwood et al. 2014, 
TEPIC - Schmidt et al. 2016) and 3) build regulatory networks to study 
a cell type of interest (Schmidt et al. 2016, Durek et al. 2016). After 
introductory presentation we will guide participants through a hands 
on practical. Therefore, we require that all participants bring their 
own laptop. Software that needs to be installed before the tutorial as 
well as data used in the tutorial will be made available on this website.

# Schedule (21.7.2017)

| Time  | Topic  |  Who |
|---|---|---|
|  2:30 - 2:45 | Introduction / gene regulation / transcription / chromatin  | Ivan G. Costa   |
|  2:45 - 3:00 | Introduction ChIP-seq peak calling  | Matthias Heinig  |
| 3:00 - 3:50  |  Practical peak calling |  Matthias Heinig |
| 4:00 - 4:15  | Coffee break  |   |
|  4:15 - 4:30 | Introduction Footprints  | Ivan G. Costa   |
|  4:30 - 4:45 |  Introduction Regulatory networks | Marcel Schulz  |
|  4:45 - 6:00 | Practical Regulatory Networks | Ivan G. Costa and Marcel Schulz |

