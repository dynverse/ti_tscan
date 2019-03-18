#!/usr/local/bin/Rscript

task <- dyncli::main()

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TSCAN)

#   ____________________________________________________________________________
#   Load data                                                               ####

counts <- task$counts
params <- task$params

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# preprocess counts
cds_prep <- TSCAN::preprocess(
  t(as.matrix(counts)),
  takelog = TRUE,
  logbase = 2,
  pseudocount = 1,
  clusternum = NULL,
  minexpr_value = params$minexpr_value,
  minexpr_percent = params$minexpr_percent,
  cvcutoff = params$cvcutoff
)

# cluster the data
cds_clus <- TSCAN::exprmclust(
  cds_prep,
  clusternum = params$clusternum,
  modelNames = params$modelNames,
  reduce = TRUE
)

# order the cells
cds_order <- TSCAN::TSCANorder(cds_clus)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# process output
pseudotime <- set_names(seq_along(cds_order), cds_order)

dimred <- cds_clus$pcareducere

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- dynwrap::wrap_data(cell_ids = rownames(dimred)) %>%
  dynwrap::add_dimred(dimred = dimred) %>%
  dynwrap::add_linear_trajectory(pseudotime = pseudotime) %>%
  dynwrap::add_timings(timings = checkpoints)

output %>% dyncli::write_output(task$output)
