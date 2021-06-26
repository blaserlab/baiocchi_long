
# set up the renv and repair with snapshot if needed
# renv::init()
# renv::snapshot()

# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/blaseRtools")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/DESeq2")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/genefilter")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/annotate")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/AnnotationDbi")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/KEGGREST")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/Biostrings")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/geneplotter")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/DoubletFinder")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/Seurat")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/SeuratDisk")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/rrvgo")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/GO.db")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/GOSemSim")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/scater")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/topGO")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/fastSave")

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

# uncomment and use the following to update the data package---------------------------------------
# latest_version <- file.info(list.files("~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/datapkg", full.names = T)) %>%
#   as_tibble(rownames = "file") %>%
#   arrange(desc(mtime)) %>%
#   dplyr::slice(1) %>%
#   pull(file) %>%
#   str_extract("datapkg/.*") %>%
#   str_replace("datapkg/", "")
#
# renv::install(paste0("~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/datapkg/",latest_version))

# use this to load the data package-------------------------------------
lazyData::requireData("baiocchi.long.datapkg")
