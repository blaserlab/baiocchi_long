
# baiocchi_long

This is the analysis project for prmt5 inhibitor single cell rna seq data.

It is currently a private repo, just like the prior one.  We will clean it up and strip out the history prior to making it public when the paper is submitted so don't worry about revealing anything.

This project has been restructured to make it easier for others (besides the creator, me) to use.  The new structure also protects the source data from unintentional changes, enhances reproducibility, and reduces system requirements.

The data now lives in a separate repository at this url:  https://github.com/blaserlab/baiocchi.long.datapkg.  But this is not where you will be getting the data from for arcane technical reasons (more below).  You will install both this analysis project and the data package onto your system.  The analysis project you can clone directly from github.  It is lightweight and you should know how to do that.

The data package is too large to be handled by github or any official repository.  So for now the working copies exist as compressed binary files (.tar.gz) on my X drive.  Once the paper is submitted we will upload the final version of the data to Dryad so reviewers can download and fiddle with it.  It will be private for peer review until the paper is accepted and then it will be public.

The best way currently to install the data package is to run the code on the file "R/dependencies.R".  

This script contains a few things worth noting.

* Commands to establish an R environment for the project.  This is like a copy of all of the packages you will be using.  If you just operated directly from your main package library, you might update something and your code stops working.  This way you keep the packages the same unless you explicitly install or update to new versions.  These are the renv:: commands

* Renv fails to automatically pull in some packages for unclear reasons, so you may have to uncomment those and manually install once.

* There is a section for loading packages used frequently in the analysis project using the standard library("package") command.

* Then there is a section that you should uncomment and run to install the latest version of the data package.  It has to be handled differently because it is not your standard package.  If you uncomment that whole block and just run it it should work.

* Finally there is a line to load the data package into your analysis project.  This is confusing because R has size limits on packages for which lazy loading is permitted.  You want lazy loading because it means that the digital objects exist only as "pointers" to an area on the hard drive until the code you run requires that object.  Only then is it loaded into memory.  It also reduces clutter in your global environment.  In order to simulate lazy loading, we use a package called lazyData.  This should be installed when you run renv::init() at the outset.  This command will load the data package into a hidden environment and then the data objects should appear when called:

```r
lazyData::requireData("baiocchi.long.datapkg")
```

After you get these things installed, your workflow would be to 

1.  pull changes in the analysis project from git
2.  check if you need to update your data package by going to the X drive and making sure your version is up to date.  If not, then uncomment and run that section in dependencies.R
3.  Source R/dependencies.R
4.  Source R/configs.R
5.  Start working on R/scratch.R.  If you have something we want to add to the data package we can. For making plots, tables and doing quick calculations we can leave those in the analysis project.

We will go over this in person too.  I know it is somewhat complex.  Good luck!
