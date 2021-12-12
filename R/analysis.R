#' ---
#' title: "Data Analysis Report"
#' author: "Brad Blaser"
#' date: "Dec 12 2021"
#' output: pdf_document
#' ---
#'
#+ include=FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo=FALSE, comment = NA)
source("R/dependencies.R")
source("R/configs.R")

#' ## Introduction and Overview
#'
#' * Data input:  ~/network/X/Labs/Blaser/single_cell/baiocchi_long_20211130
#' * Data input:  ~/network/X/Labs/Blaser/single_cell/baiocchi_long_20211209
#' * Network output:  ~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript
#'
#' Samples were processed at the CCC Genomics Shared Resource (GSR) and sequenced at Nationwide Childrens Hospital.
#'
#' Data were pre-processed using the 10X Genomics Cloud service running CellRanger v6.1.2, using GEX reference GRCh38_and_mm10-2020-A and VDJ reference vdj_GRCh38_alts_ensembl-3.1.0-3.1.0.  There were no major errors in pre-processing.
#'
#' Sample metadata provided by Mackenzie is as follows:
#'
pander::pander(analysis_configs %>% select(-pipestance_path), caption = "Sample Metadata")
#'
#' Sequencing was suboptimal for a couple of the samples due to overloading the chip (too many recovered cells).  But I don't think it is worth getting more sequencing.
#'
pander::pander(summarized_sequencing_metrics, caption = "Sequencing Metrics")

#' I used the standard functions from our blaseRtools package for qc and to remove doublets.  As before there was a population of cells with very low counts/size factor that passed initial qc that I removed with a second round of filtering.  I can provide the plots for all of these if desired and the code can be found in ```system.file("extdata/load_process.R", package = "baiocchi.long.datapackage")```.  Please note that I ran the filtering with size factor this time and renamed the cds object cds_human_pass_sf.
#'

#' ## Analysis
#' As before I applied batch correction to reduce variability in umap coordinates based on collection date.  This does not affect the underlying expression data but helps a little bit with clustering.
#'
# fix the order of the drug response  and treatment levels
colData(cds_human_pass_sf)$drug_response <- factor(colData(cds_human_pass_sf)$drug_response, levels = c("Naive", "Sensitive", "Resistant"))
colData(cds_human_pass_sf)$treatment <- factor(colData(cds_human_pass_sf)$treatment, levels = c("Untreated", "PRMT5i"))

#' The first thing I had to do was sort out all of the variables in your experimental design.  Because it is unbalanced we need to use a multivariable design in order to try to draw some conclusions, or at least make accurate hypotheses.
#'
#'  Thinking ahead to how we will analyze these data, our major response variable will be some normalized number of cells in some cell clusters which we will then define based on global transcriptional state (xy spatial distribution in UMAP space).  The major predictor we want to look at is drug response.  But since we have a complicated experimental design we need to look for known sample characteristics that might confound this interpretation and add those into the model. So we will facet the next few plots according to drug_response and other variables we need to account for, coloring cells according to sample.  Drug response will always be laid out horizontally.
#+ dev = "png", fig.align = "center", fig.width = 7.5, fig.height = 5.5
bb_var_umap(cds_human_pass_sf, "sample",
            facet_by = c("passage","drug_response"),
            rows = vars(passage),
            cols = vars(drug_response),
            foreground_alpha = 0.4,
            palette = experimental_group_palette,
            plot_title = "Drug Response x Passage")

#' The plot above facets by drug_response and passage. For the Naive and Resistant groups it looks like passage number doesn't matter in terms of xy spatial distribution relative to drug response.  All of your sensitive cells are from P6 so there is no confounding to worry about there.  Generally this looks fine in terms of drug_response and passage, however if we look at how the colors are distributed in the individual facets of the plot, it suggests unconsidered confounding variables.  For example, in the Sensitive/P6 facet you can see xy-spatial differences between green and brown.  Likewise in Resistant/P6 you can see differences between orange and blue. From what I can see in the sample names at least, we should consider some interactions between experiment and drug response and between tissue and drug response.
#'
#'  Here is the plot faceting by experiment and drug response and experiment.
#+ dev = "png", fig.align = "center", fig.width = 7.5, fig.height = 5.5
bb_var_umap(cds_human_pass_sf,
            "sample",
            facet_by = c("experiment","drug_response"),
            rows = vars(experiment),
            cols = vars(drug_response),
            foreground_alpha = 0.4,
            palette = experimental_group_palette,
            plot_title = "Drug Response x Experiment")
#'
#' Again it looks like experiment type doesn't affect the xy spatial distribution in Naive cells. However for the sensitive cells, there is a difference between sentinel and survival so we will have to consider experiment in our statistical model.


#' This plot facets by drug response and tissue.
#+ dev = "png", fig.align = "center", fig.width = 7.5, fig.height = 3.75
bb_var_umap(cds_human_pass_sf,
            "sample",
            facet_by = c("tissue","drug_response"),
            rows = vars(tissue),
            cols = vars(drug_response),
            foreground_alpha = 0.4,
            palette = experimental_group_palette,
            plot_title = "Drug Response x Tissue")

#' Here we see pretty large differences in xy spatial distribution between BM and spleen.  The BM samples look particularly problematic since (thinking ahead here) for the most part there is a pretty strong indication in the spleen cells that there are two cell populations.  The bottom one is seen in Naive and Sensitive samples and the top is seen in resistant samples.  But the problem is that the top population is pretty strong in the BM sensitive samples.  So the problem you have is that it is hard to say if the top-vs-bottom distinction more characteristic of drug response (the main question you want to answer with the experiment, so interesting) or tissue (I guess you ran BM samples because it was all you had for those conditions, not interesting).
#'
#' At this point I think it is worth developing a statistical model so we can figure out what we have and what we are going to do.  There are other variables in your sample metadata table, but these mostly align with what we have already looked at.
#'
#' First we want to fractionate the whole dataset into two groups of cells, let's put them under the heading "state" and categorize them as "resistant_like" and "sensitive_like".  Caveat being we don't yet know the true degree of association of these states with the actual samples until we do the modeling.
#'
#+ dev = "png", fig.align = "center", fig.width = 5.5, fig.height = 4.5
bb_var_umap(cds_human_pass_sf, "leiden", overwrite_labels = T)

#' Now bin them.
colData(cds_human_pass_sf)$state <- recode(colData(cds_human_pass_sf)$leiden,
                                           "11" = "sensitive_like",
                                           "10" = "sensitive_like",
                                           "15" = "sensitive_like",
                                           "16" = "sensitive_like",
                                           "2" = "sensitive_like",
                                           "6" = "sensitive_like",
                                           "8" = "sensitive_like",
                                           "12" = "sensitive_like",
                                           "4" = "sensitive_like",
                                           "1" = "resistant_like",
                                           "3" = "resistant_like",
                                           "5" = "resistant_like",
                                           "9" = "resistant_like",
                                           "7" = "resistant_like",
                                           "14" = "resistant_like",
                                           "13" = "resistant_like"
                                           )
#+ dev = "png", fig.align = "center", fig.width = 5.5, fig.height = 4.5
bb_var_umap(cds_human_pass_sf, "state", overwrite_labels = T)

#' Now for each sample calculate the log-transformed ratio betwen sensitive-like and resistant-like cells.  There is no need to control for number of cells aquired when we do it this way.  Then tack this onto the sample metadata we are interested in modeling.
#'
#' In this table, positive log2ratio indicates more resistant-like cells in the sample:

cell_state_distribution <- bb_cellmeta(cds_human_pass_sf) %>%
  group_by(sample, drug_response, tissue, experiment, state) %>%
  summarise(n = n()) %>%
  pivot_wider(everything(), values_from = n, names_from = state) %>%
  mutate(log2ratio = log2(resistant_like/sensitive_like))

pander::pander(cell_state_distribution, caption = "Cell State Distribution")


model <- lm(log2ratio ~ drug_response * experiment, data = cell_state_distribution)
# qqnorm(model$residuals)
# qqline(model$residuals)
# plot(model$residuals, model$fitted.values)

#' I fit a few different models.  The model with the best fit has experiment type as a confounding and interacting term with drug_response.
summary(model)

#' This model takes into account the different times the mice were sacrificed and whether or not they were classified as naive, sensitive or resistant.  Tissue source was not picked up in the model I think because you only had BM from the Survival samples.
#'
#' However tissue is still causing a problem.  If you look at the Estimate for drug_responseSensitive it is 3.3.  That means that after accounting for experiment type, compared to Naive, the Sensitive samples have a strong increase in the resistant-like cells, which is contrary to what we had postulated above.  It is even higher than the resistant cells.  If you look at the plot above titled Drug Response x Experiment, the brown cells are drug response Sensitive and have mostly cells in the resistant-like cluster.  These are BM cells and they basically throw off the whole analysis.
#'
#' So what we can do for starters is take out the BM cells and remodel the data.  This helps the problem.  I'm not going to give you the numbers on that to prevent confusion, because what I should to is reanalyze the whole dataset without the BM samples.
cell_state_distribution_spleen <- cell_state_distribution %>% filter(tissue == "Spleen")

spleen_model <- lm(log2ratio ~ drug_response + experiment, data = cell_state_distribution_spleen)

#' ## Conclusion
#'
#' More tbd.  It will take me a few days to redo this.  Maybe you could remind me why you used BM samples when the rest were spleen and what we lose by excluding them.  Anyways they look much different from the spleen samples so not sure there is any way to include them.
