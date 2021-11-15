################################################################################
## Adam Akullian and Allen Roberts
## Generate summary output and plots across runs within EMOD experiments
################################################################################

rm(list = ls())

library(tidyverse)
library(foreach)
library(doParallel)

## Load analysis and plotting functions
source(file.path("R", "make_results_summary_VMMC.R"))

## Set up parallels
numCores <- detectCores()
registerDoParallel(numCores)

## Suite name. This example assumes that all the runs are located within the
## same suite, but it ideally would allow arbitrary specification of 
## simulations (eg, across suites, or only certain runs within a given suite)
suite_name <- "Eswatini_Age1519_CampaignDuration_Sweep"

## Identify scenarios for a given suite
## In this test example, they are stored locally, but 
## they could be located on COMPS
scenarios <- data.frame("name" = list.files(file.path("Data", suite_name)))

## Parallel for-loop for processing and storing results in a directory specified
## by the "output_dir" argument, which is also local for this example. 
## Ideally, this output_dir could be located on a server, since these intermediate
## outputs can still be quite large.
foreach (i=1:nrow(scenarios)) %dopar% {
  
  ## Generate results summary
  make_results_summary(scenario = scenarios$name[i],
                       suite_name = suite_name,
                       data_dir = "Data",
                       output_dir = "SummaryOutput",
                       report_interval = 0.5,
                       make_plots = T,
                       make_vaccine_cohort_plots=F)
}

