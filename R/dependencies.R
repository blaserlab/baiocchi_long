
# set up the renv and repair with snapshot if needed
# renv::init()
# renv::snapshot()

# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("/usr/lib/R/site-library/blaseRtools")
# renv::install("/usr/lib/R/site-library/DESeq2")
# renv::install("/usr/lib/R/site-library/genefilter")
# renv::install("/usr/lib/R/site-library/annotate")
# renv::install("/usr/lib/R/site-library/AnnotationDbi")
# renv::install("/usr/lib/R/site-library/KEGGREST")
# renv::install("/usr/lib/R/site-library/Biostrings")
# renv::install("/usr/lib/R/site-library/geneplotter")
# renv::install("/usr/lib/R/site-library/DoubletFinder")
# renv::install("/usr/lib/R/site-library/Seurat")
# renv::install("/usr/lib/R/site-library/SeuratDisk")
# renv::install("/usr/lib/R/site-library/rrvgo")
# renv::install("/usr/lib/R/site-library/GO.db")
# renv::install("/usr/lib/R/site-library/GOSemSim")
# renv::install("/usr/lib/R/site-library/scater")
# renv::install("/usr/lib/R/site-library/topGO")
# renv::install("/usr/lib/R/site-library/fastSave")
# renv::install("ggpmisc")

# load core packages for the analysis
library("blaseRtools")
library("tidyverse")
library("monocle3")
library("circlize")
library("ComplexHeatmap")
library("lazyData")
library("cowplot")
library("RColorBrewer")
library("ggrepel")
library("ggpubr")
library("rstatix")

# uncomment and use the following to install or update the data package---------------------------------------
# bb_renv_datapkg(path = "~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/datapkg")


# use this to load the data package-------------------------------------
lazyData::requireData("baiocchi.long.datapkg")

