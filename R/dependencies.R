
# set up the renv and repair with snapshot if needed
# renv::init()

# renv::snapshot()

# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("/usr/lib/R/site-library/blaseRtools")
# renv::install("conflicted")

# load core packages for the analysis
suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("blaseRdata"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("lazyData"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
suppressPackageStartupMessages(library("conflicted"))

bb_renv_datapkg(path = "~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/datapkg")


# use this to load the data package-------------------------------------
lazyData::requireData("baiocchi.long.datapkg")

