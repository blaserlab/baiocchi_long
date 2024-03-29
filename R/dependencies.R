# package installation ----------------------------------------------------

# # By default, the newest version of all packages available  in the cache
# # at the time of project initiation are linked to the projet library

# # Use this to update the entire project library
# # to the newest versions available in the cache
# blaseRtemplates::get_new_library(newest_or_file = "newest")

# # Use this to update the entire project library
# # to another version of the project library
# blaseRtemplates::get_new_library(newest_or_file = "<path/to/file>")

# # Use this to get or update a package from the cache
# blaseRtemplates::install_one_package("blaseRtools", how = "link_from_cache")

# # If you need a new package or an update from a repository, try this:
# blaseRtemplates::install_one_package("blaseRtools", how = "new_or_update")

# # use "bioc::<package name>" for bioconductor packages
# # use "<repo/package name>" for github source packages

# # If you are installing from a "tarball", use this:
# blaseRtemplates::install_one_package("/path/to/tarball.tar.gz")

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

# load, install, and/or update the project data -----------------------------
# this requires that you are using a blaseRtemplates installation, see https://blaserlab.github.io/blaseRtemplates/articles/establish.html
# you would have to edit the path variable to point to the directory holding baiocchi.long.datapkg tar.gz file
blaseRtemplates::project_data(path = "~/network/X/Labs/Blaser/share/collaborators/baiocchi_long_manuscript/datapkg")

# alternatively, use the data function to load the data objects after installing baiocchi.long.datapkg
# courtesy of Rich Scriven https://stackoverflow.com/questions/27709936/get-a-list-of-the-data-sets-in-a-particular-package
# d <- data(package = "baiocchi.long.datapkg")
# data(list = d$results[, "Item"], package = "baiocchi.long.datapkg")
# rm(d)
