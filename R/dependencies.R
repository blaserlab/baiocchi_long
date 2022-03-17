# renv --------------------------------------------------------------------

# set up the renv from scratch

# renv::init(bioconductor = TRUE)

# restore the renv from the lockfile

# renv::restore()



# package installation ----------------------------------------------------

# # Try this first...it's faster:
# blaseRtemplates::easy_install("<package name>", how = "link_from_cache")

# # If you need a new package or an update, try this:
# blaseRtemplates::easy_install("<package name>", how = "new_or_update")

# # If you are installing from a "tarball", use this:
# blaseRtemplates::easy_install("/path/to/tarball.tar.gz")

# # use "bioc::<package name>" for bioconductor packages
# # use "<repo/package name>" for github source packages

# load core packages for the analysis -------------------------------------
suppressPackageStartupMessages(library("conflicted"))
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

bb_renv_datapkg(path = "~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/datapkg")


# use this to load the data package-------------------------------------
lazyData::requireData("baiocchi.long.datapkg")

