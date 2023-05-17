
# Resistance to PRMT5 Targeted Therapy in Mantle Cell Lymphoma

This is an analysis project.  Use this together with the companion data package to reproduce selected figures (scRNA-seq) from the manuscript.

## Steps to Reproduce

1. System Requirements
  - R v4.3
  - Rstudio
  - This software has been tested on Linux Ubuntu 18.04.6 and Windows 10
  - Loading the complete dataset occupies approximately 3.5 GB memory.

2.  Installation
  - go to the private link provided in the manuscript
  - download baiocchi.long.datapkg_<version>.tar.gz to a convenient location on your computer.
  - install the data package by clicking on the Packages tab in Rstudio, then install.  Then select install from:  package archive file and navigate/select the downloaded package.
  - clone this analysis project to your computer using git clone https://github.com/blaserlab/baiocchi_long.git
  - open the R project by double-clicking on the baiocchi_long.Rproj file
  - a list of the packages required for the project can be found in library_catalogs/blas02_baiocchi_long.tsv.  Filter for packages with status == "active".   Install these packages.
  - install custom packages from our R Universe repository using these commands:
    -  install.packages('blaseRtools', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))
    -  install.packages('blaseRtemplates', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))
    -  install.packages('blaseRdata', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))
  - edit R/dependencies.R
    - unless you are using a blaseRtemplates installation, comment out the last active line starting blaseRtemplates::project_data....
    - uncomment/activate the last 3 lines to load the data into your workspace
  - edit R/configs.R 
    - the file paths defining the figs_out  and tables_out variables should be customized for your system
  - typical time required for the first installation and data loading is approximately 10 minutes. 

3.  Instructions for use after installing and configuring
  - source R/dependencies.R
  - source R/configs.R
  - source R/manuscript_figures.R. This will generate all computationally-derived figures in the manuscript.
  - If properly configured, these scripts should run to completion in 1-2 minutes.

4.  Each computationally-generated figure panel is associated with processed data and code for visualization.  Each processed data object has its own help manual and associated processing code within the data package.  To access these resources do the following:
  - find the variable name for the panel you wish to review in R/manuscript_figures.R.
  - find the original data object used to generate that panel in the code
  - type ?<data object name> to get the help manual
  - to review processing code, go to the installed location of baiocchi.long.datapkg on your system, enter the data-raw directory and run grep --include=\*.R -rnw '.' -e "<data object name>"



