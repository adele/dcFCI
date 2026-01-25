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
  # plan("multicore", workers = n_cores) # forking
   plan("cluster", workers = n_cores)
}

scary_dataset_folder <- "../SCARY_results/dense/"
scary_output_folder <- "../SCARY_results/dense/"


net_type <- "dense"
sizes =  c("small", "med", "large") #, "x_large")
seeds = c(1,2,3)

#size = "small"
#seed = 2



run_sims = TRUE # set to FALSE if files are already processed.
restore_files = TRUE

runTrueMetrics = FALSE
runOrdMECFaithfDegree <- TRUE

# For mixed data types
# Available in R
runFCI = TRUE # FALSE # Fast Causal Inference (FCI) - 2008
runCFCI = FALSE # FALSE # Conservative FCI (cFCI) - 2012, 2014 (?)
runMIICSS = TRUE # FALSE #  MIIC_search&score (MIICSS) 2025
runDCFCI = TRUE # FALSE # Data-Compatible FCI (dcFCI) -- our work

# For only continuous (gaussian) variables
# Available in R
runBCCD = TRUE # Bayesian Constraint-Based Causal Discovery (BCCD) - 2012

# Requires pre-computation -- available in other languages
runDCD = FALSE # Differentiable Causal Discovery (DCD) - 2021
runMAGSL = FALSE # MAG Structure Learning - 2021
runGPS = FALSE  # Greedy PAG Search (GPS) - 2022

true_bic <- data.frame()
true_pag_metrics <- data.frame()
true_rmec_metrics <- data.frame()
fci_metrics <- data.frame()
cfci_metrics <- data.frame()
bccd_metrics <- data.frame()
miicss_metrics <- data.frame()
dcfci_metrics <- data.frame()
dcfci_metrics_min <- data.frame()
dcd_metrics <- data.frame()
magsl_metrics <- data.frame()
gps_metrics <- data.frame()

last_size_i = 1
last_seed_i = 0
# loading pre-computed metrics
for (size_i in 1:length(sizes)) {
  for (seed_i in 1:length(seeds)) {
    size <- sizes[size_i]
    seed <- seeds[seed_i]
    fileid <- paste0("dense_", size, "_mixed_unfaithful_", seed)
    cur_output_folder <- paste0(scary_output_folder, fileid, "/")
    cur_lat_dirs <- list.dirs(cur_output_folder, recursive = F)
    if (length(cur_lat_dirs) > 0) {
      cur_dat_dir <- paste0(cur_lat_dirs[1], "/")
      metrics_files <- list.files(cur_dat_dir,
                                  pattern=glob2rx("20260112_*.RData"),
                                  full.names = T)
      for (metric_file in metrics_files) {
        load(metric_file) # the last files have all of them together
        last_size_i = size_i
        last_seed_i = seed_i
      }
    }
  }
}


seed_i = next_seed_i = if (last_seed_i == 3) 1 else last_seed_i + 1
size_i = next_size_i = if (last_seed_i == 3) { last_size_i + 1 } else last_size_i

for (size_i in next_size_i:length(sizes)) {
  for (seed_i in next_seed_i:length(seeds)) {
    size <- sizes[size_i]
    seed <- seeds[seed_i]

    alpha = 0.05 # 0.01

    cat("\nProcessing PAG: ", size, "/", seed, "\n")

    fileid <- paste0("dense_", size, "_mixed_unfaithful_", seed)

    cur_output_folder <- paste0(scary_output_folder, fileid, "/")

    dat <- suffStat <- true_objs <- NULL
    cur_lat_dirs <- list.dirs(cur_output_folder, recursive = F)
    if (length(cur_lat_dirs) > 0) {
      cur_dat_dir <- cur_lat_dirs[1]

      dat <- read.csv(paste0(cur_dat_dir, "/", "obs_data.csv"), header = T)
      suffStat_file <- paste0(cur_dat_dir, "/", "nnGCM_suffStat.RData")
      load(suffStat_file)
      suffStat$eff_size = 0.05

      orig_suffStat_file <- paste0(cur_dat_dir, "/", "nnGCM_suffStat_orig.RData")
      if (!file.exists(orig_suffStat_file)) {
        save(suffStat, file=orig_suffStat_file)
      }

      true_objs_file <- paste0(cur_dat_dir, "/", "true_objs.RData")
      load(true_objs_file)
    }

    if (is.null(dat) || is.null(suffStat) || is.null(true_objs)) {
      stop("ERROR: SCARY dataset ", fileid, " has not been processed yet!\n")
    }

    ####################
    # Showing True PAG #
    ####################

    debug = FALSE
    true.amat.pag <- true_objs$true.amat.pag[colnames(dat), colnames(dat)]
    true.sepset <- true_objs$true.sepset
    if (debug) {
      renderAG(true.amat.pag)
      formatSepset(true.sepset)
      getSepsetResults(suffStat$citestResults, true.sepset)
    }



    ##################################
    # Assessing the true PAG metrics #
    ##################################

    indepTest <- mixedCITest
    output_folder <- paste0(cur_dat_dir, "/")

    # if (runTrueMetrics) {
    #   cur_true_bic <- NA
    #   cur_true_bic <- getBIC(true.amat.pag, dat, type="BIC")$bic
    #   true_bic <-  rbind(true_bic,
    #                        cbind(fileid=fileid, bic=cur_true_bic))
    #
    #
    #   faithfDegree <- getFaithfulnessDegree(true.amat.pag, citestResults = suffStat$citestResults)
    #
    #   true_sf_score <- getStraightforwardPAGScore2(true.amat.pag, suffStat) #  comparable, as it includes all tests
    #   true_mec_score <- getMECTargetedPAGScore(true.amat.pag, suffStat, ord=m.max,
    #                                            allowNewTests = TRUE, indepTest = mixedCITest,
    #                                            verbose = TRUE)
    #   suffStat$citestResults <- true_mec_score$citestResults
    #
    #   true_pag_metrics <- rbind(true_pag_metrics,
    #                             cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
    #                                              bic=cur_true_bic,
    #                                              sf_frechetLB=true_sf_score$frechetLB,
    #                                              sf_frechetUB=true_sf_score$frechetUB,
    #                                              sf_1mse=1 - true_sf_score$mse,
    #                                              mec_frechetLB=true_mec_score$frechetLB,
    #                                              mec_frechetUB=true_mec_score$frechetUB,
    #                                              mec_1mse=1 - true_mec_score$mse))
    # }

    if (runOrdMECFaithfDegree) {
      # for each r=0,.., p-2, computes number of tests (mec_dep and mec_indep) that are faithful
      # under alpha = 0.01 (counts1) and 0.05 (counts2), plus number of tests (counts3)
      # although not necessary, a sufficient faithfulness for dcFCI is not to have a prob greater than 50%
      # in each of the tests required in each of the order.

      rpags_mec_scores_file <- paste0(output_folder, "rpags_mec_scores_", fileid, ".RData")

      if (!file.exists(rpags_mec_scores_file)) {
        max_ord <- 2 # ncol(true.amat.pag) - 2

        selected_metrics <- c("frechetLB", "frechetUB", "mse")

        processRPAG <- function(r) {
          cat("Processing RPAG ", r, "\n")
          rpag <- getRPAG(true.amat.pag, r, true.sepset)$r.amat.pag
          rmec_score_out <- getMECTargetedPAGScore(rpag, suffStat, ord=r,
                                                   allowNewTests = TRUE,
                                                   indepTest = mixedCITest,
                                                   verbose = TRUE)
          #suffStat$citestResults <- rmec_score_out$citestResults
          r_out <- unlist(rmec_score_out)[selected_metrics]
          names(r_out) <- paste0(r, "_", names(r_out))
          return(list(r_out, citestResults=rmec_score_out$citestResults))
        }


        rpags_mec_scores <- future_lapply(0:max_ord, processRPAG,
                                           future.seed=TRUE)


        all_citestResults <- c()
        for (i in 1:length(rpags_mec_scores)) {
          all_citestResults <- rbind(all_citestResults, rpags_mec_scores[[i]]$citestResults)
        }
        citestResults <- all_citestResults[!duplicated(all_citestResults[,1:4]),]
        suffStat$citestResults <- citestResults


        save(suffStat, file=paste0(cur_dat_dir, "/", "nnGCM_suffStat.RData"))
        save(rpags_mec_scores, file=rpags_mec_scores_file)
      } else {
        load(rpags_mec_scores_file)
      }

      rpags_mec_scores_clean <- list()
      for (i in 1:length(rpags_mec_scores)) {
        rpags_mec_scores_clean[[i]] <- rpags_mec_scores[[i]][[1]]
      }

      true_rmec_metrics <- rbind(true_rmec_metrics,
                                 cbind.data.frame(fileid=fileid,
                                                  t(unlist(rpags_mec_scores_clean))))
    }

    ###############################
    # Setting up Causal Discovery #
    ###############################


    labels <- colnames(true.amat.pag)
    vars_names <- labels
    covs_names = c()
    indepTest <- mixedCITest
    fciSuffStat <- suffStat
    fciSuffStat$retall <- FALSE


    #################
    # Running dcFCI #
    #################

    if (runDCFCI) {
      #################
      # Running dcFCI #
      #################

      sel_top = 1
      #sel_top_list <- c(1,2) #,3)
      prob_sel_top = FALSE
      pH0ThreshMin = 0.2
      pH0ThreshMax = 0.9
      eff_size = 0.05

      #for (sel_top in sel_top_list) {
        dcfci_output_folder <- paste0(output_folder, "dcFCI/", "top", sel_top, "/")

        if (!file.exists(dcfci_output_folder)) {
          dir.create(dcfci_output_folder, recursive = TRUE)
        }

        eff_size_str <- if (is.null(eff_size)) "NULL" else eff_size

        dcfci_out_file <- paste0(dcfci_output_folder, "dcfci_out_", fileid, "_", eff_size_str, ".RData")

        if (!restore_files || !file.exists(dcfci_out_file)) {
          # m.max = 2; fixedGaps = NULL; fixedEdges = NULL;
          # verbose = 2; sel_top = 1; prob_sel_top = FALSE; run_parallel = TRUE;
          # allowNewTests=TRUE;
          # pH0ThreshMin=pH0ThreshMin; pH0ThreshMax=pH0ThreshMax; list.max = 500;
          # log_folder = file.path(getwd(), "tmp", "logs")
          # sapply(list.files("./R", full.names = T), source)
          done = FALSE
          while (!done) {
            cat("Processing dcFCI with pH0ThreshMax =", pH0ThreshMax, "\n")
            start_time <- proc.time()
            dcfci_out <- dcFCI(suffStat, indepTest, labels, alpha,
                               m.max = 2,
                               verbose = 2,
                               sel_top = sel_top,
                               prob_sel_top = FALSE,
                               run_parallel = TRUE,
                               allowNewTests=TRUE,
                               list.max = 1500,
                               pH0ThreshMin = pH0ThreshMin,
                               pH0ThreshMax = pH0ThreshMax,
                               order_by_mse = FALSE,
                               log_folder = file.path(dcfci_output_folder, "tmp", "logs"))
            end_time <- proc.time()
            elapsed_time <- end_time - start_time
            dcfci_out$elapsed_time <- elapsed_time

            if (!is.null(dcfci_out) && !is.null(dcfci_out$mec_score_df)) {
              save(dcfci_out, file=dcfci_out_file)
              done = TRUE
            } else {
              cat("ERROR: could not process dcFCI. Change pH0ThreshMin and/or pH0ThreshMax.")
              pH0ThreshMax = pH0ThreshMax - 0.05
            }
          }

        } else {
          load(dcfci_out_file)
        }

        ###########################################

        m.max <- dcfci_out$order_processed
        true_rmax <- getRPAG(true.amat.pag, r=m.max)
        true_rmax_pag <- true_rmax$r.amat.pag
        true_rmax_sepset <- true_rmax$r.sepset
        if (debug) {
          renderAG(true_rmax_pag)
          formatSepset(true_rmax_sepset)
          getSepsetResults(suffStat$citestResults, true_rmax_sepset)
        }
        ###########################################


        #r_true_pag <- getRPAG(true.amat.pag, r = dcfci_out$order_processed)

        dcfci_metrics_out <- getDCFCIMetrics(dcfci_out, dat,
                                       suffStat$citestResults, true_rmax_pag,
                                       checkViolations = FALSE)
        cat("\n dcFCI SHD -- mean: ", dcfci_metrics_out$dcfci_metrics_mean$shd, "min: ",
            dcfci_metrics_out$dcfci_metrics_min$shd, "\n")
        cat("\n dcFCI mec_score_up -- mean: ", dcfci_metrics_out$dcfci_metrics_mean$mec_score.2, "min: ",
            dcfci_metrics_out$dcfci_metrics_min$mec_score.2, "\n")
        if (runOrdMECFaithfDegree) {
          cat("\n true mec_score_up: ", rpags_mec_scores_clean[[m.max+1]][2], "\n")
        }


        dcfci_metrics <- rbind(dcfci_metrics,
                               cbind(fileid=fileid,
                                     eff_size=eff_size_str,
                                     sel_top = sel_top, prob_sel_top = prob_sel_top,
                                     m.max = m.max, error_message = dcfci_out$error_message,
                                     pH0ThreshMin=pH0ThreshMin,
                                     pH0ThreshMax=pH0ThreshMax,
                                     data.frame(dcfci_metrics_out$dcfci_metrics_mean)))

        dcfci_metrics_min <- rbind(dcfci_metrics_min,
                                   cbind(fileid=fileid,
                                         eff_size=eff_size_str,
                                         sel_top = sel_top, prob_sel_top = prob_sel_top,
                                         min_id = which.min(dcfci_metrics_out$dcfci_metrics$shd),
                                         m.max = m.max, error_message = dcfci_out$error_message,
                                         pH0ThreshMin=pH0ThreshMin,
                                         pH0ThreshMax=pH0ThreshMax,
                                         data.frame(dcfci_metrics_out$dcfci_metrics_min)))
      }
    #}


    #####################
    # Running baselines #
    #####################

    if (runFCI) {
      ########################
      # Running original FCI #
      ########################

      fci_output_folder <- paste0(output_folder, "FCI/")
      if (!file.exists(fci_output_folder)) {
        dir.create(fci_output_folder, recursive = TRUE)
      }

      fci_out_file <- paste0(fci_output_folder, "fci_out_", fileid, ".RData")
      if (!restore_files || !file.exists(fci_out_file)) {
        fci_out <- runFCIHelper(indepTest, fciSuffStat, alpha=alpha, m.max = m.max,
                                labels=labels, fileid=fileid,
                                output_folder=fci_output_folder,
                                savePlots = FALSE,
                                conservative=FALSE)
      } else {
        load(fci_out_file)
      }


      cur_fci_metrics <- getMetrics(true_rmax_pag, fci_out$pag,
                                    est.sepset=getPAGImpliedSepset(fci_out$pag), #fixSepsetList(fci_out$sepset),
                                    dat=dat, conservative = FALSE)


      fci_sf_score <- getStraightforwardPAGScore2(fci_out$pag, suffStat) #  comparable, as it includes all tests

      fci_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
      if (cur_fci_metrics$valid) {
        fci_mec_score <- getMECTargetedPAGScore(fci_out$pag, suffStat, ord=m.max,
                                                 allowNewTests = TRUE, indepTest = mixedCITest,
                                                 verbose = TRUE)
        suffStat$citestResults <- fci_mec_score$citestResults
      }

      fci_metrics <- rbind(fci_metrics,
                           cbind.data.frame(fileid=fileid,
                                            cur_fci_metrics,
                                            time_taken = as.numeric(fci_out$time_taken),
                                            sf_frechetLB=fci_sf_score$frechetLB,
                                            sf_frechetUB=fci_sf_score$frechetUB,
                                            sf_1mse=1 - fci_sf_score$mse,
                                            mec_frechetLB=fci_mec_score$frechetLB,
                                            mec_frechetUB=fci_mec_score$frechetUB,
                                            mec_1mse=1 - fci_mec_score$mse))
    }

    if (runCFCI) {
      ############################
      # Running conservative FCI #
      ############################

      cfci_output_folder <- paste0(output_folder, "cFCI/")
      if (!file.exists(cfci_output_folder)) {
        dir.create(cfci_output_folder, recursive = TRUE)
      }

      cfci_out_file <- paste0(cfci_output_folder, "cfci_out_", fileid, ".RData")
      if (!restore_files || !file.exists(cfci_out_file)) {
        cfci_out <- runFCIHelper(indepTest, fciSuffStat, alpha=alpha, m.max = m.max,
                                 labels=labels, fileid=fileid,
                                 output_folder=cfci_output_folder,
                                 savePlots = FALSE,
                                 conservative=TRUE)
      } else {
        load(cfci_out_file)
      }

      cur_cfci_metrics <- getMetrics(true_rmax_pag, cfci_out$pag,
                                     est.sepset=getPAGImpliedSepset(cfci_out$pag), #cfci_out$sepset,
                                     dat=dat,
                                     conservative = FALSE) # we want to know whether the PAG
      # is valid even with the "extra" circles

      cfci_sf_score <- getStraightforwardPAGScore2(cfci_out$pag, suffStat) #  comparable, as it includes all tests

      cfci_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
      if (cur_cfci_metrics$valid) {
         cfci_mec_score <- getMECTargetedPAGScore(cfci_out$pag, suffStat, ord=m.max,
                                                allowNewTests = TRUE, indepTest = mixedCITest,
                                                verbose = TRUE)
     }

      cfci_metrics <- rbind(cfci_metrics,
                            cbind.data.frame(fileid=fileid,
                                             cur_cfci_metrics,
                                             time_taken = as.numeric(cfci_out$time_taken),
                                             sf_frechetLB=cfci_sf_score$frechetLB,
                                             sf_frechetUB=cfci_sf_score$frechetUB,
                                             sf_1mse=1 - cfci_sf_score$mse,
                                             mec_frechetLB=cfci_mec_score$frechetLB,
                                             mec_frechetUB=cfci_mec_score$frechetUB,
                                             mec_1mse=1 - cfci_mec_score$mse))
    }

    if (runMIICSS) {
      # devtools::install_github("miicTeam/miicsearchscore")
      library(miicsearchscore)

      miicss_output_folder <- paste0(output_folder, "MIICS/")
      if (!file.exists(miicss_output_folder)) {
        dir.create(miicss_output_folder, recursive = TRUE)
      }

      miicss_out_file <- paste0(miicss_output_folder, "miicss_out_", fileid, ".RData")

      if (!restore_files || !file.exists(miicss_out_file)) {
        # returns the adj matrix of an AG, so we must convert it to PAG afterwards.
        miicss_out <- list()

        start_time <- Sys.time()
        adj <- tryCatch({
          run_miic_searchscore(data = dat, n_threads = n_cores)
        },
        error=function(cond) {
          print(cond)
          return(NULL)
        },
        warning=function(cond) {
          print(cond)
          return(NULL)
        })
        end_time <- Sys.time()
        miicss_out$time_taken <- end_time - start_time
        miicss_out$ag <- NULL
        miicss_out$isValid <- FALSE
        miicss_out$pag <- NULL
        if (!is.null(adj)) {
          # replace parents (2) with tails (3 in pcalg notation)
          adj[which(adj == 2, arr.ind = T)] <- 3
          # replace bidirected edges (6) with arrowheads (2 in pcalg notation)
          adj[which(adj == 6, arr.ind = T)] <- 2
          # replace children (-2) with arrowheads (2 in pcalg notation)
          adj[which(adj == -2, arr.ind = T)] <- 2

          miicss_out$isValid <- isValidAG(adj)
          if (miicss_out$isValid) {
            magg <- pcalg::pcalg2dagitty(adj, type="mag", label=colnames(adj))
            miicss_out$pag <- getTruePAG(magg)@amat
          }
        }
        save(miicss_out, file=miicss_out_file)
      } else {
        load(miicss_out_file)
      }

      #renderAG(miicss_out$pag)

      miicss_sepset <- NULL
      if (!is.null(miicss_out$pag)) {
        miicss_sepset <- getPAGImpliedSepset(miicss_out$pag)
      }

      # p <- ncol(dat)
      # complete_pag <- matrix(TRUE, nrow = p, ncol = p)
      # diag(complete_pag) <- FALSE
      # colnames(complete_pag) <- rownames(complete_pag) <- colnames(dat)
      miicss_rmax_pag <- getRPAG(miicss_out$pag, r = m.max, sepset = miicss_sepset)
      #renderAG(miicss_rmax_pag$r.amat.pag)
      #formatSepset(miicss_rmax_pag$r.sepset)

      cur_miicss_metrics <- getMetrics(true_rmax_pag, miicss_rmax_pag$r.amat.pag,
                                       est.sepset = miicss_rmax_pag$r.sepset,
                                       dat=dat,
                                       conservative = FALSE)

      miicss_sf_score <- list(frechetLB=0, frechetUB=0, mse=1)
      if (!is.null(miicss_out$pag)) {
        miicss_sf_score <- getStraightforwardPAGScore2(miicss_rmax_pag$r.amat.pag, suffStat) #  comparable, as it includes all tests
      }

      miicss_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
      if (cur_miicss_metrics$valid) {
        miicss_mec_score <- getMECTargetedPAGScore(miicss_rmax_pag$r.amat.pag, suffStat, ord=m.max,
                                                 allowNewTests = TRUE, indepTest = mixedCITest,
                                                 verbose = TRUE)
        suffStat$citestResults <- miicss_mec_score$citestResults
      }

      miicss_metrics <- rbind(miicss_metrics,
                              cbind.data.frame(fileid=fileid,
                                               cur_miicss_metrics,
                                               time_taken = as.numeric(miicss_out$time_taken),
                                               sf_frechetLB=miicss_sf_score$frechetLB,
                                               sf_frechetUB=miicss_sf_score$frechetUB,
                                               sf_1mse=1 - miicss_sf_score$mse,
                                               mec_frechetLB=miicss_mec_score$frechetLB,
                                               mec_frechetUB=miicss_mec_score$frechetUB,
                                               mec_1mse=1 - miicss_mec_score$mse))

    }

    if (runBCCD) {
      ################
      # Running BCCD #
      ################

      # devtools::install_git('https://gitlab.science.ru.nl/gbucur/RUcausal')
      library(RUcausal)

      bccd_output_folder <- paste0(output_folder, "BCCD/")
      if (!file.exists(bccd_output_folder)) {
        dir.create(bccd_output_folder, recursive = TRUE)
      }

      bccd_out_file <- paste0(bccd_output_folder, "bccd_out_", fileid, ".RData")
      if (!restore_files || !file.exists(bccd_out_file)) {
        start_time <- Sys.time()
        R <- cor(dat)
        bccd_out <- BCCD(R, nrow(dat))
        bccd_out$PAG[which(bccd_out$PAG == 3, arr.ind = T)] <- 4
        bccd_out$PAG[which(bccd_out$PAG == 1, arr.ind = T)] <- 3
        bccd_out$PAG[which(bccd_out$PAG == 4, arr.ind = T)] <- 1
        bccd_out$PAG <- t(bccd_out$PAG)
        colnames(bccd_out$PAG) <- rownames(bccd_out$PAG) <- colnames(dat)
        end_time <- Sys.time()
        bccd_out$time_taken <- end_time - start_time

        bccd_out_file <- paste0(bccd_output_folder, "bccd_out_", fileid, ".RData")
        save(bccd_out, file=bccd_out_file)
      } else {
        load(bccd_out_file)
      }

      bccd_rmax_pag <- getRPAG(bccd_out$PAG, r = m.max, sepset = getPAGImpliedSepset(bccd_out$PAG))
      #renderAG(bccd_rmax_pag$r.amat.pag)
      #formatSepset(bccd_rmax_pag$r.sepset)

      # renderAG(bccd_out$PAG, bccd_output_folder, fileid = fileid,
      #         add_index = FALSE)

      cur_bccd_metrics <- getMetrics(true_rmax_pag, bccd_rmax_pag$r.amat.pag,
                                     est.sepset = bccd_rmax_pag$r.sepset,
                                     dat=dat,
                                     conservative = FALSE) # we want to know whether the PAG
      # is valid even with the "extra" circles

      bccd_sf_score <- getStraightforwardPAGScore2(bccd_rmax_pag$r.amat.pag, suffStat) #  comparable, as it includes all tests

      bccd_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
      if (cur_bccd_metrics$valid) {
        bccd_mec_score <- getMECTargetedPAGScore(bccd_rmax_pag$r.amat.pag, suffStat, ord=m.max,
                                                   allowNewTests = TRUE, indepTest = mixedCITest,
                                                   verbose = TRUE)
        suffStat$citestResults <- bccd_mec_score$citestResults
      }

      bccd_metrics <- rbind(bccd_metrics,
                            cbind.data.frame(fileid=fileid,
                                             cur_bccd_metrics,
                                             time_taken = as.numeric(bccd_out$time_taken),
                                             sf_frechetLB=bccd_sf_score$frechetLB,
                                             sf_frechetUB=bccd_sf_score$frechetUB,
                                             sf_1mse=1 - bccd_sf_score$mse,
                                             mec_frechetLB=bccd_mec_score$frechetLB,
                                             mec_frechetUB=bccd_mec_score$frechetLB,
                                             mec_1mse=1 - bccd_mec_score$mse))
    }

    #TODO: update later to use true_rmax_pag
    if (runDCD) {
      dcd_output_file <- paste0(output_folder, "DCD2/dcd_output2.txt")
      dcd_output_time_file <- paste0(output_folder, "DCD2/dcd_output2_time.txt")
      if (file.exists(dcd_output_file) && file.exists(dcd_output_time_file)) {
        dcd_time_taken <- as.numeric(readLines(dcd_output_time_file, n = 1))

        dcd_output <- readLines(dcd_output_file, n = 3)
        di_edges <- fromJSON(dcd_output[[1]])
        bi_edges <- fromJSON(dcd_output[[2]])

        dcd_amat <- matrix(0, nrow=length(labels), ncol=length(labels))
        colnames(dcd_amat) <- rownames(dcd_amat) <- labels
        if (length(di_edges) > 0) {
          for (rowi in 1:nrow(di_edges)) {
            dcd_amat[di_edges[rowi, 2], di_edges[rowi, 1]] <- 3
            dcd_amat[di_edges[rowi, 1], di_edges[rowi, 2]] <- 2
          }
        }
        if (length(bi_edges) > 0) {
          for (rowi in 1:nrow(bi_edges)) {
            dcd_amat[bi_edges[rowi, 2], bi_edges[rowi, 1]] <- 2
            dcd_amat[bi_edges[rowi, 1], bi_edges[rowi, 2]] <- 2
          }
        }

        dcd_amag <- pcalg::pcalg2dagitty(dcd_amat, type = "mag", labels = labels)
        dcd_pag <- getTruePAG(dcd_amag)@amat
        #renderAG(dcd_pag)

        dcd_bic <- as.numeric(dcd_output[[3]])

        cur_dcd_metrics <- getMetrics(true.amat.pag, dcd_pag,
                                      est.sepset = getPAGImpliedSepset(dcd_pag),
                                      dat=dat,
                                      conservative = FALSE)

        dcd_sf_score <- getStraightforwardPAGScore2(dcd_pag, suffStat) #  comparable, as it includes all tests

        dcd_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
        if (cur_dcd_metrics$valid) {
          dcd_mec_score <- getMECTargetedPAGScore(dcd_pag, suffStat, ord=Inf)
        }

        dcd_metrics <- rbind(dcd_metrics,
                             cbind.data.frame(fileid=fileid,
                                              cur_dcd_metrics,
                                              dcd_bic = dcd_bic,
                                              time_taken = dcd_time_taken,
                                              sf_frechetLB=dcd_sf_score$frechetLB,
                                              sf_frechetUB=dcd_sf_score$frechetUB,
                                              sf_1mse=1 - dcd_sf_score$mse,
                                              mec_frechetLB=dcd_mec_score$frechetLB,
                                              mec_frechetUB=dcd_mec_score$frechetUB,
                                              mec_1mse=1 - dcd_mec_score$mse))
      } else {
        stop(paste0("Error: dcd output file not found: ", dcd_output_file))
      }
    }


    #TODO: update later to use true_rmax_pag
    if (runMAGSL) {
      library(readr)
      magsl_output_file <- paste0(output_folder, "MAGSL2/magsl_output_bb2.txt")
      magsl_scorer_time_file <- paste0(output_folder, "MAGSL2/scorer_time.txt")
      magsl_search_time_file <- paste0(output_folder, "MAGSL2/search_time.txt")
      if (file.exists(magsl_output_file) && file.exists(magsl_scorer_time_file) &&
          file.exists(magsl_search_time_file)) {
        magsl_output <- read_tsv(magsl_output_file, col_names = FALSE)
        magsl_amat.mag <- as.matrix(magsl_output[1:length(labels),1:length(labels)])
        magsl_out_score <- as.numeric(magsl_output[length(labels)+1, 1])
        colnames(magsl_amat.mag) <- rownames(magsl_amat.mag) <- labels

        magsl_amat.mag <- pcalg::pcalg2dagitty(magsl_amat.mag, type = "mag", labels = labels)
        magsl_amat.pag <- getTruePAG(magsl_amat.mag)@amat
        #renderAG(magsl_amat.pag)

        magsl_time_taken <- as.numeric(readLines(magsl_scorer_time_file, n = 1)) +
          as.numeric(readLines(magsl_search_time_file, n = 1))


        cur_magsl_metrics <- getMetrics(true.amat.pag, magsl_amat.pag,
                                        est.sepset = getPAGImpliedSepset(magsl_amat.pag),
                                        dat=dat,
                                        conservative = FALSE)


        magsl_sf_score <- getStraightforwardPAGScore2(magsl_amat.pag, suffStat) #  comparable, as it includes all tests

        magsl_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
        if (cur_magsl_metrics$valid) {
          magsl_mec_score <- getMECTargetedPAGScore(magsl_amat.pag, suffStat, ord=Inf)
        }

        magsl_metrics <- rbind(magsl_metrics,
                               cbind.data.frame(fileid=fileid,
                                                cur_magsl_metrics,
                                                magsl_score = magsl_out_score,
                                                time_taken = magsl_time_taken,
                                                sf_frechetLB=magsl_sf_score$frechetLB,
                                                sf_frechetUB=magsl_sf_score$frechetUB,
                                                sf_1mse=1 - magsl_sf_score$mse,
                                                mec_frechetLB=magsl_mec_score$frechetLB,
                                                mec_frechetUB=magsl_mec_score$frechetUB,
                                                mec_1mse=1 - magsl_mec_score$mse))
      } else {
        stop(paste0("Error: magsl output file not found: ", magsl_output_file))
      }
    }

    #TODO: update later to use true_rmax_pag
    if (runGPS) {
      gps_output_file <- paste0(output_folder_sim, "GPS2/gps_output.txt")
      if (file.exists(gps_output_file)) {
        gps_output <- read.csv(gps_output_file, header = FALSE)
        gps_amat.pag <- as.matrix(gps_output[1:(nrow(gps_output)-2), ])
        gps_out_score <- gps_output[nrow(gps_output)-1, 1]
        gps_time_taken <- gps_output[nrow(gps_output), 1]
        colnames(gps_amat.pag) <- rownames(gps_amat.pag) <- labels
        #renderAG(gps_amat.pag)

        cur_gps_metrics <- getMetrics(true.amat.pag, gps_amat.pag,
                                      est.sepset = getPAGImpliedSepset(gps_amat.pag),
                                      dat=dat,
                                      conservative = FALSE)

        gps_sf_score <- getStraightforwardPAGScore2(gps_amat.pag, suffStat) #  comparable, as it includes all tests

        gps_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
        if (cur_gps_metrics$valid) {
          gps_mec_score <- getMECTargetedPAGScore(gps_amat.pag, suffStat, ord=Inf)
        }

        gps_metrics <- rbind(gps_metrics,
                             cbind.data.frame(fileid=fileid,
                                              cur_gps_metrics,
                                              gps_score = gps_out_score,
                                              time_taken = gps_time_taken,
                                              sf_frechetLB=gps_sf_score$frechetLB,
                                              sf_frechetUB=gps_sf_score$frechetUB,
                                              sf_1mse=1 - gps_sf_score$mse,
                                              mec_frechetLB=gps_mec_score$frechetLB,
                                              mec_frechetUB=gps_mec_score$frechetUB,
                                              mec_1mse=1 - gps_mec_score$mse))
      } else {
        stop(paste0("Error: gps output file not found: ", gps_output_file))
      }
    }

    #######################
    # Saving result files #
    #######################

    save(suffStat, file=paste0(cur_dat_dir, "/", "nnGCM_suffStat.RData"))

    if (runOrdMECFaithfDegree) {
      save(true_rmec_metrics, file = paste0(output_folder, "20260112_true_rmec_metrics.RData"))
    }

    # if (runTrueMetrics) {
    #   save(true_bic, file = paste0(output_folder, "20260112_true_bic.RData"))
    #   save(true_pag_metrics, file = paste0(output_folder, "20260112_true_pag_metrics.RData"))
    # }

    if (runFCI) {
      save(fci_metrics, file = paste0(output_folder, "20260112_fci_metrics.RData"))
    }

    if (runCFCI) {
      save(cfci_metrics, file = paste0(output_folder, "20260112_cfci_metrics.RData"))
    }

    if (runMIICSS) {
      save(miicss_metrics, file = paste0(output_folder, "20260112_miicss_metrics.RData"))
    }

    if (runBCCD) {
      save(bccd_metrics, file = paste0(output_folder, "20260112_bccd_metrics.RData"))
    }

    if (runDCD) {
      save(dcd_metrics, file = paste0(output_folder, "20260112_dcd_metrics.RData"))
    }

    if (runMAGSL) {
      save(magsl_metrics, file = paste0(output_folder, "20260112_magsl_metrics.RData"))
    }

    if (runGPS) {
      save(gps_metrics, file = paste0(output_folder, "20260112_gps_metrics.RData"))
    }

    if (runDCFCI) {
      save(dcfci_metrics,
           file = paste0(output_folder, "20260112_dcfci_metrics.RData"))
      save(dcfci_metrics_min,
           file = paste0(output_folder, "20260112_dcfci_metrics_min.RData"))
    }

  } # end seeds loop
} # end sizes loop

