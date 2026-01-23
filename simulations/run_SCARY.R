rm(list=ls())

library(FCI.Utils)
library(dcFCI)
library(jsonlite)

library(doFuture)
library(future.apply)
library(rje)
library(dplyr)
source("./simulations/sim_utils.R")
source("./simulations/helper_functions.R")
source("./R/metrics_PAG.R")

# If running in parallel:
run_parallel = TRUE

if (run_parallel) {
  require(doFuture)
  require(future.apply)
  n_cores <- parallel::detectCores() # or change to the number of cores required
  # plan("multisession", workers = n_cores) # for running on RStudio
  # plan("multisession", workers = n_cores)
  # plan("multicore", workers = n_cores) # forking
  plan("cluster", workers = n_cores)
}

type = "dense"
#sizes = c("small", "med", 
sizes = c("large") #, "x_large")
seeds = c(2,3)


data_folder <- "../../SCARY_CLeaR2023/SCARY-main/CSV_unfaithful/sparse/"
if (type == "dense") {
  data_folder <- "../../SCARY_CLeaR2023/SCARY-main/CSV_unfaithful/dense/"
}

size = "small"
seed = 2
type = "dense"

output_folder <- "../../SCARY_results/sparse/"
if (type == "dense") {
  output_folder <- "../../SCARY_results/dense/"
}

for (size in sizes) {
  for (seed in seeds) {
    cat("Processing size:", size, ", seed:", seed, "\n")

    data_id <- paste0(size, "_mixed_unfaithful_", seed)
    if (type == "dense") {
      data_id <- paste0("dense_", size, "_mixed_unfaithful_", seed)
    }

    cur_output_folder <- paste0(output_folder, data_id, "/")
    if (!dir.exists(cur_output_folder)) {
      dir.create(cur_output_folder, recursive = T)
    }

    cur_lat_dirs <- list.dirs(cur_output_folder, recursive = F)
    if (length(cur_lat_dirs) > 0) {
      cur_dat_dir <- cur_lat_dirs[1]

      suffStat_file <- paste0(cur_dat_dir, "/", "nnGCM_suffStat.RData")
      if (!file.exists(suffStat_file)) {
        data <- read.csv(paste0(cur_dat_dir, "/", "obs_data.csv"), header = T)
        dim(data)

        indepTest <- mixedCITest
        vars_names <- colnames(data) # obs variables / nodes of the graph
        suffStat <- getMixedCISuffStat(dat = data, vars_names = vars_names,
                                       covs_names = c(),
                                       method = "nnGCM", verbose = TRUE)
        suffStat$retall = TRUE  # Required if computing probs for dcFCI
        suffStat$nruns = 5
        suffStat$compute_MC_pvalue = FALSE

        vars_df <- data[,vars_names, drop=FALSE]

        # Compute all conditional independence tests in parallel
        citestResults <- getAllCITestResults(
          vars_df, indepTest, suffStat, m.max = 2,
          computeProbs = TRUE,        # Required if running dcFCI
          eff_size = 0.05,            # Minimum effect size threshold when computeProbs = TRUE
          fileid = data_id,           # Unique identifier for the dataset / analysis
          recover_citestResults = TRUE, # Set TRUE to resume computation after a crash
          results_folder = cur_dat_dir,
          log_folder = cur_dat_dir
        )

        #ret <- indepTest(3, 8, NULL, suffStat = suffStat)
        suffStat$citestResults <- citestResults
        save(suffStat, file=suffStat_file)
      } else {
        cat("suffStat has been already processed.\n")
      }
    } else {
      cat("ERROR: Data folder not available!\n")
    }
  }
}






