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
  n_cores <- 12
  # plan("multisession", workers = n_cores)
  plan("multicore", workers = n_cores) # forking
  #plan("cluster", workers = n_cores)
}


# ############
# # To check #
# ############
#
# data_type = "mixed" #  # option in {continuous, mixed}
# N <- 200
# pag_id <- 1
# sim <- 10

# dcfci_metrics[which(dcfci_metrics$shd > 0 &
#   dcfci_metrics$mec_score.2 == true_pag_metrics$mec_frechetUB),]

#########################
# Simulation Parameters #
#########################

sim_ids = 1:30
sample_sizes = c(200, 500) #, 1000, 5000, 10000, 50000)


sim_true_pags_file <- paste0("./simulations/random_p5_true_pags_list.RData")
load(sim_true_pags_file)
#lapply(true_pags_list, renderAG)

pag_ids <- 1:length(true_pags_list)

restore_files = TRUE
restore_suffStat = TRUE
debug = FALSE

####################


run_sims = TRUE # set to FALSE if files are already processed.

# For mixed data types
# Available in R
runFCI = FALSE # Fast Causal Inference (FCI) - 2008
runCFCI = FALSE # Conservative FCI (cFCI) - 2012, 2014 (?)
runMIICSS = FALSE #  MIIC_search&score (MIICSS) 2025
runDCFCI = TRUE # Data-Compatible FCI (dcFCI) -- our work

# For only continuous (gaussian) variables
# Available in R
runBCCD = FALSE # Bayesian Constraint-Based Causal Discovery (BCCD) - 2012
# Requires pre-computation -- available in other languages
runDCD = FALSE # Differentiable Causal Discovery (DCD) - 2021
runMAGSL = FALSE # MAG Structure Learning - 2021
runGPS = FALSE  # Greedy PAG Search (GPS) - 2022


runTrueMetrics = TRUE
runOrdMECFaithfDegree <- TRUE

run_plots = FALSE

# check_processed: TRUE, skips directly;
#                  FALSE, loads pre-computed files and computes metrics again
check_processed = FALSE

if (run_sims) {
  #######################
  # Running Simulations #
  #######################

  for (data_type in c("mixed", "continuous")) {
    output_folder <- paste0("../dcFCI_Simulations/",
                            data_type, "/")

    if (data_type == "continuous") {
      runTrueMetrics = FALSE
      runOrdMECFaithfDegree <- FALSE
    } else {
      runTrueMetrics = TRUE
      runOrdMECFaithfDegree <- TRUE
    }

    n_processed = 0
    if (check_processed) {
      metrics_files <- list.files(output_folder,
                                  pattern=glob2rx("20260112_*.RData"))
      for (metric_file in metrics_files) {
        load(paste0(output_folder, metric_file))
      }
      # Note: since all files are saved in the end, they all have the same
      # number of rows in case of crash.
      n_processed = dim(subset(dcfci_metrics, sel_top == 2))[1]
      #n_processed = dim(fci_metrics)[1]
    } else {
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
    }

    res_i = 0
    for (N in sample_sizes) {
      for (pag_id in pag_ids) {
        for (sim in sim_ids) {
          res_i = res_i + 1
          if (res_i <= n_processed) {
            next
          }

          alpha = 0.05 # 0.01

          cat("\n\nPAG: ", pag_id, "N: ", N, "sim: ", sim, "\n")

          fileid <- paste0(c(data_type, pag_id, N, sim), collapse="_")


          ####################
          # Showing True PAG #
          ####################

          true.amat.pag <- true_pags_list[[pag_id]]
          true.sepset <- getPAGImpliedSepset(true.amat.pag)
          if (debug) {
           renderAG(true.amat.pag)
           formatSepset(true.sepset)
          }

          ###################
          # Loading dataset #
          ###################

          dat_folder <- paste0("../dcFCI_Simulations/",
                               data_type, "/",  sprintf("%02d", pag_id),
                               "/n", N, "/sim_", sim, "/")

          dat_file <- list.files(dat_folder, pattern=glob2rx("dat_*.csv"))
          if (length(dat_file) == 0) {
            message("Data not found for ", fileid)
            next
          }


          dat <- read.csv(file=paste0(dat_folder, dat_file))
          cat("\nDataset with", ncol(dat), "columns.\n")

          cat_cols <- as.numeric(which(sapply(dat, is.character)))
          dat[, cat_cols] <- lapply(dat[, cat_cols], as.factor)


          #if (runMAGSL) {
            magsl_csv <- paste0(dat_folder, "/", "magsl_", dat_file)
            if (!file.exists(magsl_csv)) {
              cat("Saving dataset compatible with MAGSL...\n")
              write.table(dat, file=magsl_csv,
                          row.names = FALSE, col.names = FALSE, sep = " ")
            }
          #}


          eff_size = if (data_type == "continuous") 0.05 else 0.05 # suffStat (1)  was 0.01 for mixed
          eff_size_str <- if (is.null(eff_size)) "NULL" else eff_size
          if (is.null(eff_size)) {
            suffStat_file <- paste0(dat_folder, "suffStat_sim_", sim, "_", eff_size_str, ".RData")
          } else { #eff_zize = 0.05
            suffStat_file <- paste0(dat_folder, "suffStat_sim_", sim, ".RData")
            if (data_type == "mixed") { # && N >= 1000) {
              # for N=200 and N=500:
              #   suffStat2 is with eff_size = 0.05
              #   suffStat uses eff_size = 0.01
              # for N >= 1000:
              #   suffStat2 is with eff_size = 0.01
              #   suffStat uses eff_size = 0.05
              suffStat_file <- paste0(dat_folder, "suffStat2_sim_", sim, ".RData")
            }
          }

          if (restore_suffStat && file.exists(suffStat_file)) {
            load(suffStat_file)
          } else {
            #dat$D <- factor(dat$D)
            #dat$E <- factor(dat$E)

            vars_names <- colnames(true.amat.pag) # only observed variables
            covs_names = c()
            indepTest <- mixedCITest

            suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
            vars_df <- dat[,vars_names, drop=FALSE]
            citestResults <- getAllCITestResults(vars_df, indepTest, suffStat,
                                                 m.max=Inf, computeProbs = TRUE,
                                                 eff_size = eff_size)
            suffStat$citestResults <- citestResults
            save(suffStat, file=suffStat_file)
          }

          ##################################
          # Assessing the true PAG metrics #
          ##################################

          if (runTrueMetrics) {
            cur_true_bic <- NA
            if (data_type == "continuous") {
              cur_true_bic <- getBIC(true.amat.pag, dat, type="BIC")$bic
              true_bic <-  rbind(true_bic,
                                 cbind(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
                                       bic=cur_true_bic))
            }

            true_sf_score <- getStraightforwardPAGScore2(true.amat.pag, suffStat) #  comparable, as it includes all tests
            true_mec_score <- getMECTargetedPAGScore(true.amat.pag, suffStat, ord=Inf)

            true_pag_metrics <- rbind(true_pag_metrics,
                                  cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
                                                   bic=cur_true_bic,
                                                   sf_frechetLB=true_sf_score$frechetLB,
                                                   sf_frechetUB=true_sf_score$frechetUB,
                                                   sf_1mse=1 - true_sf_score$mse,
                                                   mec_frechetLB=true_mec_score$frechetLB,
                                                   mec_frechetUB=true_mec_score$frechetUB,
                                                   mec_1mse=1 - true_mec_score$mse))
          }

          if (runOrdMECFaithfDegree) {
            # for each r=0,.., p-2, computes number of tests (mec_dep and mec_indep) that are faithful
            # under alpha = 0.01 (counts1) and 0.05 (counts2), plus number of tests (counts3)
            # although not necessary, a sufficient faithfulness for dcFCI is not to have a prob greater than 50%
            # in each of the tests required in each of the order.
            max_ord <- ncol(true.amat.pag) - 2

            selected_metrics <- c("frechetLB", "frechetUB", "mse")
            rpags_mec_scores <- lapply(0:max_ord, function(r) {
              rpag <- getRPAG(true.amat.pag, r, true.sepset)$r.amat.pag
              r_out <- unlist(getMECTargetedPAGScore(rpag, suffStat, ord=r))[selected_metrics]
              names(r_out) <- paste0(r, "_", names(r_out))
              r_out
            })

            true_rmec_metrics <- rbind(true_rmec_metrics,
                                        cbind.data.frame(data_type=data_type, N=N,
                                                         pag_id=pag_id, sim=sim,
                                                         t(unlist(rpags_mec_scores))))
          }


          ###############################
          # Setting up Causal Discovery #
          ###############################

          output_folder <- paste0("../dcFCI_Simulations/",
                                  data_type, "/")
          output_folder_pagid <- paste0(output_folder,  sprintf("%02d", pag_id), "/")
          output_folder_N <- paste0(output_folder_pagid, "n", N, "/")
          output_folder_sim <- paste0(output_folder_N, "sim_", sim, "/")


          if (output_folder_sim != dat_folder) {
            #head(dat)
            write.csv(dat, file=paste0(output_folder_sim, dat_file),
                      row.names = FALSE)
            save(suffStat, file=paste0(output_folder_sim, "suffStat_sim_", sim, ".RData"))
          }


          labels <- colnames(true.amat.pag)
          vars_names <- labels
          covs_names = c()
          indepTest <- mixedCITest

          #####################
          # Running baselines #
          #####################

          if (runFCI) {
            ########################
            # Running original FCI #
            ########################

            fci_output_folder <- paste0(output_folder_sim, "FCI2/")
            if (!file.exists(fci_output_folder)) {
              dir.create(fci_output_folder, recursive = TRUE)
            }

            fci_out_file <- paste0(fci_output_folder, "fci_out_", fileid, ".RData")
            if (!restore_files || !file.exists(fci_out_file)) {
              fci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                                      labels=labels, fileid=fileid,
                                      output_folder=fci_output_folder,
                                      savePlots = FALSE,
                                      conservative=FALSE)
            } else {
              load(fci_out_file)
            }

            eval_dat <- if (data_type == "continuous") dat else NULL
            cur_fci_metrics <- getMetrics(true.amat.pag, fci_out$pag,
                                          est.sepset=getPAGImpliedSepset(fci_out$pag), #fixSepsetList(fci_out$sepset),
                                          dat=eval_dat, conservative = FALSE)


            fci_sf_score <- getStraightforwardPAGScore2(fci_out$pag, suffStat) #  comparable, as it includes all tests

            fci_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
            if (cur_fci_metrics$valid) {
              fci_mec_score <- getMECTargetedPAGScore(fci_out$pag, suffStat, ord=Inf)
            }

            fci_metrics <- rbind(fci_metrics,
                                      cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
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

            cfci_output_folder <- paste0(output_folder_sim, "cFCI2/")
            if (!file.exists(cfci_output_folder)) {
              dir.create(cfci_output_folder, recursive = TRUE)
            }

            cfci_out_file <- paste0(cfci_output_folder, "cfci_out_", fileid, ".RData")
            if (!restore_files || !file.exists(cfci_out_file)) {
              cfci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                                       labels=labels, fileid=fileid,
                                       output_folder=cfci_output_folder,
                                       savePlots = FALSE,
                                       conservative=TRUE)
            } else {
              load(cfci_out_file)
            }

            eval_dat <- if (data_type == "continuous") dat else NULL
            cur_cfci_metrics <- getMetrics(true.amat.pag, cfci_out$pag,
                                           est.sepset=getPAGImpliedSepset(cfci_out$pag), #cfci_out$sepset,
                                           dat=eval_dat,
                                           conservative = FALSE) # we want to know whether the PAG
                                                                   # is valid even with the "extra" circles

            cfci_sf_score <- getStraightforwardPAGScore2(cfci_out$pag, suffStat) #  comparable, as it includes all tests

            cfci_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
            if (cur_cfci_metrics$valid) {
              cfci_mec_score <- getMECTargetedPAGScore(cfci_out$pag, suffStat, ord=Inf)
            }

            cfci_metrics <- rbind(cfci_metrics,
                                 cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
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

            miicss_output_folder <- paste0(output_folder_sim, "MIICSS/")
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

            eval_dat <- if (data_type == "continuous") dat else NULL
            cur_miicss_metrics <- getMetrics(true.amat.pag, miicss_out$pag,
                                           est.sepset = miicss_sepset,
                                           dat=eval_dat,
                                           conservative = FALSE)

            miicss_sf_score <- list(frechetLB=0, frechetUB=0, mse=1)
            if (!is.null(miicss_out$pag)) {
              miicss_sf_score <- getStraightforwardPAGScore2(miicss_out$pag, suffStat) #  comparable, as it includes all tests
            }

            miicss_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
            if (cur_miicss_metrics$valid) {
              miicss_mec_score <- getMECTargetedPAGScore(miicss_out$pag, suffStat, ord=Inf)
            }

            miicss_metrics <- rbind(miicss_metrics,
                                  cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
                                                   cur_miicss_metrics,
                                                   time_taken = as.numeric(miicss_out$time_taken),
                                                   sf_frechetLB=miicss_sf_score$frechetLB,
                                                   sf_frechetUB=miicss_sf_score$frechetUB,
                                                   sf_1mse=1 - miicss_sf_score$mse,
                                                   mec_frechetLB=miicss_mec_score$frechetLB,
                                                   mec_frechetUB=miicss_mec_score$frechetUB,
                                                   mec_1mse=1 - miicss_mec_score$mse))

          }

          if (runBCCD && data_type == "continuous") {
            ################
            # Running BCCD #
            ################

            # devtools::install_git('https://gitlab.science.ru.nl/gbucur/RUcausal')
            library(RUcausal)

            bccd_output_folder <- paste0(output_folder_sim, "BCCD2/")
            if (!file.exists(bccd_output_folder)) {
              dir.create(bccd_output_folder, recursive = TRUE)
            }

            bccd_out_file <- paste0(bccd_output_folder, "bccd_out_", fileid, ".RData")
            if (!restore_files || !file.exists(bccd_out_file)) {
              start_time <- Sys.time()
              R <- cor(dat)
              bccd_out <- BCCD(R, N)
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

            # renderAG(bccd_out$PAG, bccd_output_folder, fileid = fileid,
            #         add_index = FALSE)

            cur_bccd_metrics <- getMetrics(true.amat.pag, bccd_out$PAG,
                                           est.sepset = getPAGImpliedSepset(bccd_out$PAG),
                                           dat=dat,
                                           conservative = FALSE) # we want to know whether the PAG
            # is valid even with the "extra" circles

            bccd_sf_score <- getStraightforwardPAGScore2(bccd_out$PAG, suffStat) #  comparable, as it includes all tests

            bccd_mec_score <- list(frechetLB=0, frechetUB=0, mse=1)
            if (cur_bccd_metrics$valid) {
              bccd_mec_score <- getMECTargetedPAGScore(bccd_out$PAG, suffStat, ord=Inf)
            }
            bccd_metrics <- rbind(bccd_metrics,
                                  cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
                                                   cur_bccd_metrics,
                                                   time_taken = as.numeric(bccd_out$time_taken),
                                                   sf_frechetLB=bccd_sf_score$frechetLB,
                                                   sf_frechetUB=bccd_sf_score$frechetUB,
                                                   sf_1mse=1 - bccd_sf_score$mse,
                                                   mec_frechetLB=bccd_mec_score$frechetLB,
                                                   mec_frechetUB=bccd_mec_score$frechetLB,
                                                   mec_1mse=1 - bccd_mec_score$mse))

          }

          if (runDCD && data_type == "continuous") {
            dcd_output_file <- paste0(output_folder_sim, "DCD2/dcd_output2.txt")
            dcd_output_time_file <- paste0(output_folder_sim, "DCD2/dcd_output2_time.txt")
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
                                   cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
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

          if (runMAGSL && data_type == "continuous") {
            library(readr)
            magsl_output_file <- paste0(output_folder_sim, "MAGSL2/magsl_output_bb2.txt")
            magsl_scorer_time_file <- paste0(output_folder_sim, "MAGSL2/scorer_time.txt")
            magsl_search_time_file <- paste0(output_folder_sim, "MAGSL2/search_time.txt")
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
                                   cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
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

          if (runGPS && data_type == "continuous") {
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
                                     cbind.data.frame(data_type=data_type, N=N, pag_id=pag_id, sim=sim,
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


          #################
          # Running dcFCI #
          #################

          if (runDCFCI) {
            #################
            # Running dcFCI #
            #################

            sel_top_list <- c(1,2) #,3)
            prob_sel_top = FALSE
            pH0ThreshMin = 0.3

            for (sel_top in sel_top_list) {
              # dcFCI4 - see 4t results - check_structure = FALSE, with sel_top 1 & 2, prob_sel_top = 0,
              # dcFCI5 - see 5t results - using mec_score, check_structure = FALSE, with sel_top 1 & 2, prob_sel_top = 0,

              if (data_type == "continuous") {
                dcfci_output_folder <- paste0(output_folder_sim, "dcFCI6/", "top", sel_top, "/")
              } else {
                dcfci_output_folder <- paste0(output_folder_sim, "dcFCI7/", "top", sel_top, "/")
              }
              if (!file.exists(dcfci_output_folder)) {
                dir.create(dcfci_output_folder, recursive = TRUE)
              }

              eff_size_str <- if (is.null(eff_size)) "NULL" else eff_size

              dcfci_out_file <- paste0(dcfci_output_folder, "dcfci_out_", fileid, "_", eff_size_str, ".RData")
              if (data_type == "continuous") {
                dcfci_out_file_old <- paste0(dcfci_output_folder, "dcfci_out_", fileid, ".RData")
                if (!file.exists(dcfci_out_file) && file.exists(dcfci_out_file_old)) {
                  file.rename(dcfci_out_file_old, dcfci_out_file)
                }
              }

              redo_date <- as.POSIXct("2026-01-17 15:00:00", tz="CET")

              if (!restore_files || !file.exists(dcfci_out_file) ||
                  (N < 1000 && as.POSIXct(file.info(dcfci_out_file)$ctime) < redo_date)) {
                # m.max = Inf; fixedGaps = NULL; fixedEdges = NULL;
                # verbose = 2; sel_top = 1; prob_sel_top = FALSE; run_parallel = TRUE;
                # allowNewTests=TRUE; pH0ThreshMin=0.3; pH0ThreshMax=1; list.max = 500;
                # log_folder = file.path(getwd(), "tmp", "logs")

                start_time <- proc.time()
                dcfci_out <- dcFCI(suffStat, indepTest, labels, alpha,
                                   m.max = Inf,
                                   verbose = FALSE,
                                   sel_top = sel_top,
                                   prob_sel_top = FALSE,
                                   run_parallel = TRUE,
                                   allowNewTests=TRUE,
                                   list.max = 1000,
                                   pH0ThreshMin=pH0ThreshMin,
                                   pH0ThreshMax = 1,
                                   log_folder = file.path(dcfci_output_folder, "tmp", "logs"))
                end_time <- proc.time()
                elapsed_time <- end_time - start_time
                dcfci_out$elapsed_time <- elapsed_time

                if (!is.null(dcfci_out)) {
                  save(dcfci_out, file=dcfci_out_file)
                }
              } else {
                load(dcfci_out_file)

                if (!is.null(dcfci_out$fit_dcfci)) {
                  dcfci_out_out <- dcfci_out
                  dcfci_out <- dcfci_out$fit_dcfci
                  dcfci_out$elapsed_time <- dcfci_out_out$time_taken
                  save(dcfci_out, file=dcfci_out_file)
                }
              }

              eval_dat <- if (data_type == "continuous") dat else NULL
              metrics_out <- getDCFCIMetrics(dcfci_out, eval_dat,
                                             suffStat$citestResults, true.amat.pag)
              cat("\n dcFCI SHD: ", metrics_out$dcfci_metrics$shd, "\n")
              cat("\n dcFCI mec_score_up: ", metrics_out$dcfci_metrics$mec_score.2, "\n")
              if (runTrueMetrics) {
                cat("\n true mec_score_up: ", true_mec_score$frechetUB, "\n")
              }


              dcfci_metrics <- rbind(dcfci_metrics,
                                     cbind(data_type=data_type, N=N, pag_id=pag_id, sim=sim, eff_size=eff_size_str,
                                           sel_top = sel_top, prob_sel_top = prob_sel_top,
                                           pH0ThreshMin=pH0ThreshMin,
                                           data.frame(metrics_out$dcfci_metrics)))

              dcfci_metrics_min <- rbind(dcfci_metrics_min,
                                         cbind(data_type=data_type, N=N, pag_id=pag_id, sim=sim, eff_size=eff_size_str,
                                               sel_top = sel_top, prob_sel_top = prob_sel_top,
                                               pH0ThreshMin=pH0ThreshMin,
                                               data.frame(metrics_out$dcfci_metrics)))
            }
          }


          # if (runFCI && runCFCI && runDCFCI) {
          #   cat("SHD\n")
          #   print(tail(cbind.data.frame(fci=fci_metrics$shd, cfci=cfci_metrics$shd, #bccd=bccd_metrics$shd,
          #           dcfci1=dcfci_metrics[[1]]$shd, dcfci2=dcfci_metrics[[2]]$shd)))
          #   cat("FDR\n")
          #   print(tail(cbind.data.frame(fci=fci_metrics$fdr, cfci=cfci_metrics$fdr, #bccd=bccd_metrics$fdr,
          #     dcfci1=dcfci_metrics[[1]]$fdr, dcfci2=dcfci_metrics[[2]]$fdr)))
          #   cat("FOR\n")
          #   print(tail(cbind.data.frame(fci=fci_metrics$fomr, cfci=cfci_metrics$fomr, #bccd=bccd_metrics$fomr,
          #     dcfci1=dcfci_metrics[[1]]$fomr, dcfci2=dcfci_metrics[[2]]$fomr)))
          # }


          #######################
          # Saving result files #
          #######################

          if (runOrdMECFaithfDegree) {
            save(true_rmec_metrics, file = paste0(output_folder, "20260112_true_rmec_metrics.RData"))
          }

          if (runTrueMetrics) {
            if (data_type == "continuous") {
              save(true_bic, file = paste0(output_folder, "20260112_true_bic.RData"))
            }
            save(true_pag_metrics, file = paste0(output_folder, "20260112_true_pag_metrics.RData"))
          }

          if (runFCI) {
            save(fci_metrics, file = paste0(output_folder, "20260112_fci_metrics.RData"))
          }

          if (runCFCI) {
            save(cfci_metrics, file = paste0(output_folder, "20260112_cfci_metrics.RData"))
          }

          if (runMIICSS) {
            save(miicss_metrics, file = paste0(output_folder, "20260112_miicss_metrics.RData"))
          }

          if (runBCCD && data_type == "continuous") {
            save(bccd_metrics, file = paste0(output_folder, "20260112_bccd_metrics.RData"))
          }

          if (runDCD && data_type == "continuous") {
            save(dcd_metrics, file = paste0(output_folder, "20260112_dcd_metrics.RData"))
          }

          if (runMAGSL && data_type == "continuous") {
            save(magsl_metrics, file = paste0(output_folder, "20260112_magsl_metrics.RData"))
          }

          if (runGPS && data_type == "continuous") {
            save(gps_metrics, file = paste0(output_folder, "20260112_gps_metrics.RData"))
          }

          if (runDCFCI) {
              save(dcfci_metrics,
                   file = paste0(output_folder, "20260112_dcfci_metrics.RData"))
              save(dcfci_metrics_min,
                   file = paste0(output_folder, "20260112_dcfci_metrics_min.RData"))
          }
        } # end sim loop
      } # end N loop
    } # end pag_id loop
  } # end data_type loop
}



run_plots <- FALSE
if (run_plots) {
  # test for H1: mu_A > mu_B
  wilcox_summary <- function(grpA, grpB) {
    wtest <- wilcox.test(grpA, # group A
                         grpB, # group B
                         paired=TRUE,  #  paired Wilcoxon signed-rank test
                         # alternative = "less", # H1: mu_A < mu_B
                         alternative = "two.sided", #  "greater", # H0: mu_B < mu_A; H1: mu_A > mu_B
                         conf.int=FALSE, conf.level=0.95) # confidence interval for the
    return(wtest)
  }


  #data_type = "continuous"
  data_type = "mixed" #


  sel_top = 1

  # How much do I lose not checking the top 2 in terms of finding the true PAG?
  # What is the strategy? Check boths and check how the top scored PAGs change?


  output_folder <- paste0("../dcFCI_Simulations/", data_type, "/")

  metrics_files <- list.files(output_folder,
                              pattern=glob2rx("20260112_*.RData"))
  for (metric_file in metrics_files) {
    load(paste0(output_folder, metric_file))
  }

#
#   if (data_type == "continuous") {
#     load(paste0(output_folder, "true_bic3.RData"))
#     load(paste0(output_folder, "bccd_metrics3t.RData"))
#     load(paste0(output_folder, "dcd_metrics3t.RData"))
#     load(paste0(output_folder, "magsl_metrics3t.RData"))
#     load(paste0(output_folder, "gps_metrics3t.RData"))
#   }
#   load(paste0(output_folder, "faithf_degree3.RData"))
#   load(paste0(output_folder, "rmec_faithf_degree3.RData"))
#   load(paste0(output_folder, "true_metrics3.RData"))
#   load(paste0(output_folder, "fci_metrics3t.RData"))
#   load(paste0(output_folder, "cfci_metrics3t.RData"))


 for (sel_top in 1:2) {

  dcfci_metrics8 <- data.frame()
  dcfci_metrics_min8 <- data.frame()
  improved_ids <- c()
  worsened_ids <- c()


  # Loading dcfci_metrics2
  load(paste0(output_folder, "dcfci_metrics6t_sel_top_", sel_top, ".RData"))
  dcfci_metrics6 <- dcfci_metrics2
  # Loading dcfci_metrics_min2
  load(paste0(output_folder, "dcfci_metrics_min6t_sel_top_", sel_top, ".RData"))
  dcfci_metrics_min6 <- dcfci_metrics_min2

  # Loading dcfci_metrics2
  load(paste0(output_folder, "dcfci_metrics7t_sel_top_", sel_top, ".RData"))
  dcfci_metrics7 <- dcfci_metrics2
  # Loading dcfci_metrics_min2
  load(paste0(output_folder, "dcfci_metrics_min7t_sel_top_", sel_top, ".RData"))
  dcfci_metrics_min7 <- dcfci_metrics_min2

  for (i in 1:1200) {
    if (is.na(dcfci_metrics7$shd[i]) ||
        dcfci_metrics_min6$shd[i] < dcfci_metrics_min7$shd[i]) {
      if (!is.na(dcfci_metrics_min7$shd[i])) {
        worsened_ids <- c(worsened_ids, i)
      }
      dcfci_metrics8 <- rbind(dcfci_metrics8,
                                  dcfci_metrics6[i, , drop=FALSE])
      dcfci_metrics_min8 <- rbind(dcfci_metrics_min8,
                                      dcfci_metrics_min6[i, , drop=FALSE])
    } else {
      if (dcfci_metrics_min6$shd[i] > dcfci_metrics_min7$shd[i]) {
        improved_ids <- c(improved_ids, i)
        print(c(dcfci_metrics_min6$shd[i], dcfci_metrics_min7$shd[i]))
      }
      dcfci_metrics8 <- rbind(dcfci_metrics8,
                                  dcfci_metrics7[i, , drop=FALSE])
      dcfci_metrics_min8 <- rbind(dcfci_metrics_min8,
                                      dcfci_metrics_min7[i, , drop=FALSE])
    }
  }
  dcfci_metrics2 <- dcfci_metrics8
  dcfci_metrics_min2 <- dcfci_metrics_min8

  save(dcfci_metrics2, file = paste0(output_folder,
                                     "dcfci_metrics8t_sel_top_", sel_top, ".RData"))
  save(dcfci_metrics_min2, file = paste0(output_folder,
                                     "dcfci_metrics_min8t_sel_top_", sel_top, ".RData"))
 }


  sel_top = 1
  # Loading dcfci_metrics2
  load(paste0(output_folder, "dcfci_metrics8t_sel_top_", sel_top, ".RData"))
  dcfci_metrics_k1 <- dcfci_metrics2
  # Loading dcfci_metrics_min2
  load(paste0(output_folder, "dcfci_metrics_min8t_sel_top_", sel_top, ".RData"))
  dcfci_metrics_k1_min <- dcfci_metrics_min2

  sel_top = 2
  # Loading dcfci_metrics2
  load(paste0(output_folder, "dcfci_metrics8t_sel_top_", sel_top, ".RData"))
  dcfci_metrics_k2 <- dcfci_metrics2
  # Loading dcfci_metrics_min2
  load(paste0(output_folder, "dcfci_metrics_min8t_sel_top_", sel_top, ".RData"))
  dcfci_metrics_k2_min <- dcfci_metrics_min2

  sel_top = 2
  for (sel_top in c(1,2)) {

  # Loading dcfci_metrics2
  load(paste0(output_folder, "dcfci_metrics8t_sel_top_", sel_top, ".RData"))

  # Loading dcfci_metrics_min2
  load(paste0(output_folder, "dcfci_metrics_min8t_sel_top_", sel_top, ".RData"))

  dcfci_metrics <- dcfci_metrics_min2

  dcfci_metrics$pag_id <- as.numeric(dcfci_metrics$pag_id)
  dcfci_metrics$N <- as.numeric(dcfci_metrics$N)
  length_agg <- aggregate(dcfci_metrics[, "shd"],
                          by = as.list(dcfci_metrics[, c("N", "pag_id")]), FUN=length)

  # checking if the execution for all 30 sims/pag for all 10 pags are done
  good_pag_ids <- as.numeric(names(which(table(
    subset(length_agg, x == 30)$pag_id) == length(sample_sizes))))
  length(good_pag_ids)


  ########################
  # Specifying baselines #
  ########################

  methods = c("fci", "cfci")
  methods_labels <- c("FCI", "cFCI")
  if (data_type == "continuous") {
    methods <- c(methods, "bccd", "dcd",
                 "magsl",
                 "gps")
    methods_labels <- c(methods_labels, "BCCD", "DCD",
                        "MAGSL",
                        "GPS")
  }
  methods <- c(methods, "dcfci")
  methods_labels <- c(methods_labels, "dcFCI")

  metrics_labels <- c("SHD", "FDR", "FOR", "Compatibility Score", "Time Taken")
  metrics = c("shd", "fdr", "fomr", "trivial_score.2", "time_taken")
  if (data_type == "continuous") {
    metrics_labels <- c(metrics_labels, "BIC", "BIC Distance")
    metrics <- c(metrics, "bic", "bic_dist")
  }

  library(RColorBrewer)
  n <- length(methods_labels)
  method_colors = brewer.pal(n = n, name = "Set3")


  figures_folder <- paste0(output_folder, "figures_final8/")

  if (!file.exists(figures_folder)) {
    dir.create(figures_folder, recursive = TRUE)
  }


  tables_folder <- paste0(output_folder, "tables_final8/")

  if (!file.exists(tables_folder)) {
    dir.create(tables_folder, recursive = TRUE)
  }

  methods_labels_k1k2 <- c(methods_labels[1:(length(methods_labels)-1)],
                           "dcFCI_1", "dcFCI_2")


  ######################################################################
  # Show that we can obtain the true PAG even without adj faithfulness #
  ######################################################################

  # What is the shd when faithfulness = 1? And other types of faithfulness?
  # What is the faithfulness of those with shd = 0, and with fdr = 0

  # Selects all that have at least 50% of trivial_score
  # Show only true PAG recovery rate
  # table with the number faithful data according to ...
  # and then how is the prop of recovery rate in each of them.
  true_faith_df <- data.frame()
  for(cur_N in sample_sizes) {
    cur_N_ids <- which(true_metrics$N == cur_N)

    all_metric_lists <- list(fci_metrics[cur_N_ids, , drop = FALSE],
                             cfci_metrics[cur_N_ids, , drop = FALSE])

    if (data_type == "continuous") {
      all_metric_lists[[length(all_metric_lists) + 1]] <- bccd_metrics[cur_N_ids, , drop = FALSE]
      all_metric_lists[[length(all_metric_lists) + 1]] <- dcd_metrics[cur_N_ids, , drop = FALSE]
      if ("magsl" %in% methods) {
        all_metric_lists[[length(all_metric_lists) + 1]] <- magsl_metrics[cur_N_ids, , drop = FALSE]
      }
      all_metric_lists[[length(all_metric_lists) + 1]] <- gps_metrics[cur_N_ids, , drop = FALSE]
    }
    all_metric_lists[[length(all_metric_lists) + 1]] <-  dcfci_metrics_k1_min[cur_N_ids, , drop = FALSE]
    all_metric_lists[[length(all_metric_lists) + 1]] <-  dcfci_metrics_k2_min[cur_N_ids, , drop = FALSE]

    cur_faithf_degree <- subset(faithf_degree, N == cur_N)
    cur_true_metrics <- subset(true_metrics, N == cur_N)
    cur_rmec_faithf_degree <- subset(rmec_faithf_degree, N == cur_N)

    f_faithf_ids <- which(cur_faithf_degree$f_dep.2 / cur_faithf_degree$f_dep.3 == 1 &
                            cur_faithf_degree$f_indep.2 / cur_faithf_degree$f_indep.3 == 1)
    b_faithf_ids <- which(cur_true_metrics$trivial_score.2 >= 0.5)
    #a_faithf_ids <- which(cur_faithf_degree$fskel_dep.2 / cur_faithf_degree$fskel_dep.3 == 1 &
    #                      cur_faithf_degree$fskel_indep.2 / cur_faithf_degree$fskel_indep.3 == 1)
    rmec_faithf_ids <- which(cur_rmec_faithf_degree$`0_mec_score2` >= 0.5 &
                               cur_rmec_faithf_degree$`1_mec_score2` >= 0.5 &
                               cur_rmec_faithf_degree$`2_mec_score2` >= 0.5 &
                               cur_rmec_faithf_degree$`3_mec_score2` >= 0.5 )
    #mec_faithf_ids <- which(cur_rmec_faithf_degree$`0_mec_score2` >= 0.5)

    u_ids <- setdiff(1:nrow(cur_true_metrics),
                         unique(c(f_faithf_ids, b_faithf_ids, rmec_faithf_ids)))

    for (metrics_list in all_metric_lists) {
      #n_true_a_faithf <- length(which(metrics_list[a_faithf_ids, "shd"] == 0))
      #n_fdr0_a_faithf <- length(which(metrics_list[a_faithf_ids, "fdr"] == 0))
      n_true_f_faithf <- length(which(metrics_list[f_faithf_ids, "shd"] == 0))
      n_true_b_faithf <- length(which(metrics_list[b_faithf_ids, "shd"] == 0))
      n_true_rmec_faithf <- length(which(metrics_list[rmec_faithf_ids, "shd"] == 0))
      #n_true_mec_faithf <- length(which(metrics_list[mec_faithf_ids, "shd"] == 0))
      n_true_unfaithf <- length(which(metrics_list[u_ids, "shd"] == 0))


      true_faith_df <- rbind.data.frame(true_faith_df,
                                       c(N=cur_N,
                                         #n_afaith = length(a_faithf_ids),
                                         #n_atrue = n_true_a_faithf,
                                         #n_afdr0 = n_fdr0_a_faithf,
                                         n_ffaith = length(f_faithf_ids),
                                         n_ftrue = n_true_f_faithf,
                                         n_bfaith = length(b_faithf_ids),
                                         n_btrue = n_true_b_faithf,
                                         n_rmecfaith =length(rmec_faithf_ids),
                                         n_rmectrue = n_true_rmec_faithf,
                                         #n_mecfaith = length(mec_faithf_ids),
                                         #n_mectrue = n_true_mec_faithf,
                                         n_unf = length(u_ids),
                                         n_unftrue = n_true_unfaithf))
    }
  }
  true_faith_df <-  cbind(method=rep(methods_labels_k1k2, length(sample_sizes)), true_faith_df)
  colnames(true_faith_df) <- c("Algorithm", "N",
                               "n_ffaith", "ff_recov",
                               "n_bfaith", "bf_recov",
                               "n_rmecfaith", "rmecf_recov",
                               "n_unfaith", "unf_recov")
  true_faith_df

  library(xtable)
  print(xtable(true_faith_df, floating=FALSE, latex.environments=NULL,
               display=c("s", "s", rep("d", 9)),
               digits=c(0,0, rep(3,9))),
        math.style.exponents = TRUE,  include.rownames=FALSE)


  perc_true_faith_df <- true_faith_df

  perc_true_faith_df[,4] <- signif((true_faith_df[,4]/true_faith_df[,3]) * 100, 3)
  perc_true_faith_df[,6] <- signif((true_faith_df[,6]/true_faith_df[,5]) * 100, 3)
  perc_true_faith_df[,8] <- signif((true_faith_df[,8]/true_faith_df[,7]) * 100, 3)
  perc_true_faith_df[,10] <- signif((true_faith_df[,10]/true_faith_df[,9]) * 100, 3)
  perc_true_faith_df

  print(xtable(true_faith_df[, c(1,2,3,5)], floating=FALSE, latex.environments=NULL,
               display=c("s","s","d","d","d"),
               digits=c(0,0,3,3,3)),
        math.style.exponents = TRUE,  include.rownames=FALSE)



  # # Proportion of adjacency faithful datasets
  #
  # distr_adj_score <- data.frame()
  # for (cur_N in sample_sizes) {
  #   cur_true_metrics <- subset(true_metrics, N == cur_N)
  #   min_score_cat <- cut(cur_true_metrics$adj_score.2,
  #                        breaks = seq(0, 1, length.out = 11))
  #   distr_adj_score <- rbind.data.frame(distr_adj_score,
  #                                      c(table(min_score_cat)))
  # }
  # colnames(distr_adj_score) <- levels(min_score_cat)
  # distr_adj_score <- cbind(N = sample_sizes, distr_adj_score)
  # distr_adj_score <- cbind(distr_adj_score,
  #                          faithful = apply(distr_adj_score[, c(7:11)], 1, sum))
  #
  # distr_adj_score
  #
  #
  # # Proportion of mec (adj + ori) faithful datasets
  # distr_mec_score <- data.frame()
  # for (cur_N in sample_sizes) {
  #   cur_true_metrics <- subset(true_metrics, N == cur_N)
  #   min_score_cat <- cut(cur_true_metrics$mec_score.2,
  #                        breaks = seq(0, 1, length.out = 11))
  #   distr_mec_score <- rbind.data.frame(distr_mec_score,
  #                                       c(table(min_score_cat)))
  # }
  # colnames(distr_mec_score) <- levels(min_score_cat)
  # distr_mec_score <- cbind(N = sample_sizes, distr_mec_score)
  #
  # distr_mec_score <- cbind(distr_mec_score,
  #                          faithful = apply(distr_mec_score[, c(7:11)], 1, sum))
  # distr_mec_score




  ###############################################
  # Checking Number of PAGs in dcFCI's PAG List #
  ###############################################

  ntop_df <- dcfci_metrics[, c("N", "pag_id", "sim", "ntop")]
  ntop_df <- as.data.frame(lapply(ntop_df , as.numeric))
  ntop_df <- ntop_df[order(ntop_df$N), ]

  n_top_pags <- list()
  for (cur_N in sample_sizes) {
    n_top_pags[[as.character(cur_N)]] <- table(subset(ntop_df, N == cur_N)$ntop)
  }
  n_top_pags

  sink(file=paste0(tables_folder, "n_top_pags_", data_type, "_sel_top_", sel_top, ".txt"))
  print(n_top_pags)
  sink()

  #################################
  # Checking number of violations #
  #################################


  # If valid, then checks if compatible. So, we maybe we should compute number
  # of instances, not percentage...
  viol_df <- data.frame()
  for(cur_N in sample_sizes) {
    all_metric_lists <- list( subset(fci_metrics, N == cur_N & pag_id %in% good_pag_ids),
                              subset(cfci_metrics, N == cur_N & pag_id %in% good_pag_ids))
    if (data_type == "continuous") {
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(bccd_metrics, N == cur_N & pag_id %in% good_pag_ids)
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcd_metrics, N == cur_N & pag_id %in% good_pag_ids)
      if ("magsl" %in% methods) {
        all_metric_lists[[length(all_metric_lists) + 1]] <- subset(magsl_metrics, N == cur_N & pag_id %in% good_pag_ids)
      }
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(gps_metrics, N == cur_N & pag_id %in% good_pag_ids)
    }
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcfci_metrics, N == cur_N & pag_id %in% good_pag_ids)

    for (metrics_list in all_metric_lists) {
        # viol_df <- rbind.data.frame(viol_df,  c(N=cur_N, colMeans(
        #   metrics_list[, c("valid", "viol", "msepViol")], na.rm=TRUE)))
         viol_df <- rbind.data.frame(viol_df,  c(N=cur_N, colSums(
           metrics_list[, c("valid", "viol")], na.rm=TRUE)
           ))
    }
    # viol_df <- rbind(viol_df, c(N=cur_N, 1, 0, 0))
  }
  viol_df <-  cbind(method=rep(methods_labels, length(sample_sizes)), viol_df)
  colnames(viol_df) <- c("method", "N", "valid", "viol")

  viol_df[, "invalid"] <-  length(true_pags_list) * 30 - viol_df[, "valid"]
  viol_df

  viol_df2 <- viol_df[, c(1,2,3,5)]
  viol_df2 <- reshape::melt(viol_df2, id.vars = c("N", "method"))
  viol_df2 <- reshape(viol_df2, direction="wide", idvar = c("N", "variable"), timevar = c( "method"))
  viol_df2 <- viol_df2[order(viol_df2$N), ]
  viol_df2

  library(xtable)
  sink(file=paste0(tables_folder, "viol_df2_", data_type, "_sel_top_", sel_top, ".txt"))
  print(xtable(viol_df2, floating=FALSE, latex.environments=NULL,
               display=c("s", "d", "s", rep("d", length(methods))),
               digits=c(0,3,0,rep(3, length(methods)))),
        math.style.exponents = TRUE,  include.rownames=FALSE)
  sink()





  ###################################
  # Checking True PAG Recovery Rate #
  ###################################

  # Note: using dcfci_metrics, so true pag is also the only one recovered in the list
  n_true_df <- data.frame()
  for (cur_N in sample_sizes) {
    cur_faithf_degree <- subset(faithf_degree, N == cur_N)
    cur_true_metrics <- subset(true_metrics, N == cur_N)

    #plot(cur_faithf_degree$f_dep.1/cur_faithf_degree$f_dep.3, cur_true_metrics$trivial_score.2)
    #plot(cur_faithf_degree$f_dep.1/cur_faithf_degree$f_dep.3, metrics_list$shd)
    #plot(cur_true_metrics$mec_score.2, metrics_list$shd)


    all_metric_lists <- list( subset(fci_metrics, N == cur_N & pag_id %in% good_pag_ids),
                              subset(cfci_metrics, N == cur_N & pag_id %in% good_pag_ids))
    if (data_type == "continuous") {
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(bccd_metrics, N == cur_N & pag_id %in% good_pag_ids)
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcd_metrics, N == cur_N & pag_id %in% good_pag_ids)
      if ("magsl" %in% methods) {
        all_metric_lists[[length(all_metric_lists) + 1]] <- subset(magsl_metrics, N == cur_N & pag_id %in% good_pag_ids)
      }
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(gps_metrics, N == cur_N & pag_id %in% good_pag_ids)
    }
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcfci_metrics_k1, N == cur_N & pag_id %in% good_pag_ids)
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcfci_metrics_k2, N == cur_N & pag_id %in% good_pag_ids)

    for (metrics_list in all_metric_lists) {
      total_pags <- nrow(metrics_list)
      if ("truePAGInd" %in% colnames(metrics_list)) {
        id_0_PAGs = nrow(subset(metrics_list, truePAGInd == 1))
        probid_0_PAGs = nrow(subset(metrics_list, truePAGProbInd == 1))
      } else {
        id_0_PAGs = NA
        probid_0_PAGs = NA
      }
      # shd_0_ids <- which(metrics_list$shd == 0)
      # faithf_degree_shd0 <- cur_faithf_degree[shd_0_ids, ]
      # summary(faithf_degree_shd0$f_dep.1 /faithf_degree_shd0$f_dep.3)
      # summary(faithf_degree_shd0$fadj_dep.1 /faithf_degree_shd0$fadj_dep.3)
      # summary(faithf_degree_shd0$fmec_dep.1 /faithf_degree_shd0$fmec_dep.3)
      # summary(faithf_degree_shd0$fskel_dep.1 /faithf_degree_shd0$fskel_dep.3)
      # summary(cur_true_metrics[shd_0_ids, ]$adj_score.2)
      # summary(cur_true_metrics[shd_0_ids, ]$trivial_score.2)
      # summary(cur_true_metrics[shd_0_ids, ]$mec_score.2)



      n_true_df <- rbind.data.frame(n_true_df, data.frame(N=cur_N,
                                                          SHD_0_PAGs = nrow(subset(metrics_list, shd == 0)),
                                                          FDR_FOMR_0_PAGs = nrow(subset(metrics_list, fdr == 0 & fomr == 0)),
                                                          FDR_0_PAGs = nrow(subset(metrics_list, fdr == 0)),
                                                          id_0_PAGs = id_0_PAGs,
                                                          probid_0_PAGs = probid_0_PAGs))
                     #perc=round(nrow(subset(metrics_list, shd == 0))/total_pags, 2)))

    }
  }
  n_true_df <- cbind(method=rep(methods_labels_k1k2, length(sample_sizes)), n_true_df)
  n_true_df

  n_true_values_df <- n_true_df[, c(3,6,7)]
  perc_true_values_df <- t(apply(n_true_values_df, 1, function(x) {
    round(x/300 * 100, 2) }))

  n_true_values_df <- apply(n_true_values_df, 1, function (x) {
    if (any(is.na(x))) { x[1] }
    else {paste0(x[1], " (", x[2]," / ", x[3], ")")}
  })

  perc_true_values_df <- apply(perc_true_values_df, 1, function (x) {
    if (any(is.na(x))) { x[1] }
    else {paste0(x[1], " (", x[2]," / ", x[3], ")")}
  })

  n_true_perc_df <- cbind(n_true_df[, 1:2], n_true_values_df, perc_true_values_df)

  n_true_perc_df <- reshape::melt(n_true_perc_df, id.vars = c("N", "method"))
  n_true_perc_df <- reshape(n_true_perc_df, direction="wide", idvar = c("N", "variable"), timevar = c( "method"))
  n_true_perc_df <- n_true_perc_df[order(n_true_perc_df$N), ]
  n_true_perc_df

  n_true_perc_df$variable <- factor(n_true_perc_df$variable,
                                    labels = c("counts", "perc"))
  colnames(n_true_perc_df) <- c("N", "Recovered", methods_labels_k1k2)

  sink(file=paste0(tables_folder, "true_pag_rec_rate_",
                   data_type, "_sel_top_", sel_top, ".txt"))
  print(xtable(n_true_perc_df,
               floating=FALSE, latex.environments=NULL,
               display=c("s", "d", "s", rep("s", length(methods_labels_k1k2))),
               digits=c(0, 3, 0, rep(0, length(methods_labels_k1k2)))),
        math.style.exponents = TRUE,  include.rownames=FALSE)
  sink()


  #########################################################
  # Comparing dcFCI and baselines w.r.t. all metrics      #
  # only_valid: if TRUE, compares only valid PAGs. This.  #
  # comparison demonstrates that our score help           #
  # beyond the validation checks                          #
  #########################################################

  library(DescTools)

  only_valid = FALSE
  for (only_valid in c(TRUE, FALSE)) {
    metrics_long_df <- c()
    pvalues_df <- data.frame()
    for(cur_N in sample_sizes) {
      cur_dcfci_metrics <-  subset(dcfci_metrics, N == cur_N & pag_id %in% good_pag_ids)
      all_metric_lists <- list( subset(fci_metrics, N == cur_N & pag_id %in% good_pag_ids),
                                subset(cfci_metrics, N == cur_N & pag_id %in% good_pag_ids))
      if (data_type == "continuous") {
        all_metric_lists[[length(all_metric_lists) + 1]] <- subset(bccd_metrics, N == cur_N & pag_id %in% good_pag_ids)
        all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcd_metrics, N == cur_N & pag_id %in% good_pag_ids)
        if ("magsl" %in% methods) {
          all_metric_lists[[length(all_metric_lists) + 1]] <- subset(magsl_metrics, N == cur_N & pag_id %in% good_pag_ids)
        }
        all_metric_lists[[length(all_metric_lists) + 1]] <- subset(gps_metrics, N == cur_N & pag_id %in% good_pag_ids)
      }
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcfci_metrics, N == cur_N & pag_id %in% good_pag_ids)

      for (metric in metrics) {
        if (!only_valid && !metric %in% metrics[c(1:3, 5)]) {
          next
        }
        for (metrics_list_i in 1:(length(methods))) {
          cur_method <- methods[metrics_list_i]
          metrics_list <- all_metric_lists[[metrics_list_i]]
          cur_valid_compat <- 1:nrow(metrics_list)
          #cur_valid_compat <- which(cur_dcfci_metrics$ntop == 1)
          if (only_valid) {
            cur_valid_compat <- which(metrics_list$valid == TRUE & metrics_list$viol == FALSE)
          }
          len <- length(cur_valid_compat)
          grpA <- as.numeric(metrics_list[cur_valid_compat, metric])
          grpB <- as.numeric(cur_dcfci_metrics[cur_valid_compat, metric])
          if (metric == "bic_dist") {
            grpA <- as.numeric(metrics_list[cur_valid_compat, "bic"])
            grpB <- as.numeric(cur_dcfci_metrics[cur_valid_compat, "bic"])
            cur_true_bic_metrics <- subset(true_bic, N == cur_N & pag_id %in% good_pag_ids)
            cur_true_bic <- as.numeric(cur_true_bic_metrics[cur_valid_compat, "bic"])
            grpA = abs(grpA - cur_true_bic)
            grpB = abs(grpB - cur_true_bic)
          }
          #stats_A <- signif(summary(grpA)[c(2,3,5)],2)
          #stats_B <- signif(summary(grpB)[c(2,3,5)],2)
          stats_A <- signif(summary(grpA),2)
          stats_B <- signif(summary(grpB),2)
          #pv1 <- wilcox_summary(grpA, grpB)$p.value # gprA > grpB
          #pv2 <- wilcox_summary(grpB, grpA)$p.value # gprB > grpA
          if (cur_method == "dcfci") {
            pvalue = " "
          } else {
            pvalue <- signif(SignTest(x=grpA, y =grpB, conf.level = 0.95, alternative = "two.sided")$p.value,2)
          }
          decision <- "--"
          #if (pvalue < 0.05 && stats_A[2] != stats_B[2]) {
          #    decision <- if (stats_A[2] < stats_B[2]) "uparrow" else "downarrow"
          #}
          if (pvalue < 0.05) {
              if (stats_A[3] < stats_B[3] ||
                  stats_A[4] < stats_B[4]) decision = "uparrow"
              else if (stats_A[3] > stats_B[3] ||
                       stats_A[4] > stats_B[4]) decision = "downarrow"
              #else if (stats_A[2] > stats_B[2]) decision = "downarrow"
          }

          #p_sign_l <- SignTest(x=grpA, y =grpB, conf.level = 0.95, alternative = "less")
          #p_sign_g <- SignTest(x=grpA, y =grpB, conf.level = 0.95, alternative = "greater")

          pvalues_df <- rbind.data.frame(pvalues_df,
                                      data.frame(N=cur_N, metric=metric, len=len,
                                                 t(as.numeric(stats_A[1:6])),
                                                 #t(as.numeric(stats_B)),
                                                 #pvalue=pv1,
                                                 pvalue=pvalue,
                                                 decision
                                                 #p_sign_l=p_sign_l$p.value,
                                                 #p_sign_g=p_sign_g$p.value
                                      ))
          metrics_long_df <- rbind.data.frame(metrics_long_df,
                                           cbind.data.frame(val=grpA, N=cur_N,
                                                            metric=metric, method=cur_method))
          # if (metrics_list_i == length(all_metric_lists)) {
          #   cur_method <- "dcfci"
          #   metrics_long_df <- rbind.data.frame(metrics_long_df,
          #                                    cbind.data.frame(val=grpB, N=cur_N,  metric=metric, method=cur_method))
          # }
        }
      }
    }
    n_metrics <- length(metrics)
    if (!only_valid) {
      n_metrics <- 4
    }

    pvalues_df = cbind(
      method=rep(methods_labels[1:(length(methods_labels))], length(sample_sizes) * n_metrics), pvalues_df)

    pvalues_df$method <- factor(pvalues_df$method,
                                levels = methods_labels[1:(length(methods_labels))],
                                labels = methods_labels[1:(length(methods_labels))])

    if (!only_valid) {
      pvalues_df <- pvalues_df[, -which(colnames(pvalues_df) == "len")]
      colnames(pvalues_df) <- c("Algorithm", "N", "metric",
                                "Min","Q1", "Med", "Mean", "Q3", "Max",
                                "pvalue", "decision")
    } else {
      colnames(pvalues_df) <- c("Algorithm", "N", "metric", "nsims",
                                "Min","Q1", "Med", "Mean", "Q3", "Max",
                                "pvalue", "decision")
    }


    head(pvalues_df, 7)


    metrics_long_df$method <- factor(metrics_long_df$method, levels=methods,
                                  labels=methods_labels)

    colnames(metrics_long_df) <- c("val", "N", "metric", "Algorithm")


    library(ggplot2)
    for (cur_metric_i in 1:5) {
      if (!only_valid && !cur_metric %in% metrics[c(1:3, 5)]) {
        next
      }
      cur_metric <- metrics[cur_metric_i]
      cur_metrics_long <- subset(metrics_long_df, metric == cur_metric) # & Algorithm != "DCD")
      p2 <- ggplot(cur_metrics_long, aes(x=Algorithm,y=val,fill=Algorithm))+
        geom_boxplot() +
        scale_fill_manual(values=method_colors) +
        ggtitle(paste0("Boxplot of ", metrics_labels[cur_metric_i], " by Algorithm and Sample Size")) +
        ylab(metrics_labels[cur_metric_i]) + xlab("Algorithm") +
        facet_wrap(~N)
      p2
      ggsave(paste0(figures_folder, "boxplots_only_valid_", only_valid, "_", cur_metric,
                    "_", data_type, "_sel_top_", sel_top, ".png"),
             width=7, height =5)
    }

    if (only_valid && data_type == "continuous") {
      for (cur_metric_i in 6:7) {
        cur_metric <- metrics[cur_metric_i]
        cur_metrics_long <- subset(metrics_long_df, metric == cur_metric)
        p2 <- ggplot(cur_metrics_long, aes(x=Algorithm,y=val,fill=Algorithm))+
          geom_boxplot() +
          scale_fill_manual(values=method_colors) +
          ggtitle(paste0("Boxplot of ", metrics_labels[cur_metric_i], " by Algorithm and Sample Size")) +
          ylab(metrics_labels[cur_metric_i]) + xlab("Algorithm") +
          facet_wrap(~N, scales = "free")
        p2
        ggsave(paste0(figures_folder, "boxplots_only_valid_", only_valid, "_", cur_metric,
                      "_", data_type, "_sel_top_", sel_top, ".png"),
               width=7, height =5)
      }
    }


    #############################
    # Generating p-value Tables #
    #############################

    for (cur_metric in metrics) {
      if (!only_valid && !cur_metric %in% metrics[c(1:3, 5)]) {
        next
      }

      if (!only_valid) {
        display_vec <- c("s", "s", "d",  rep("g", 7), "s")
        digits_vec <- c(0, 0, rep(2, 8), 0)
      } else {
        display_vec <- c("s", "s", "d", "d", rep("g", 7), "s")
        digits_vec <- c(0, 0, rep(2, 9), 0)
      }
      cur_metric_pvalues_df <- subset(pvalues_df, metric == cur_metric)[ , -which(colnames(pvalues_df) == "metric")]
      cur_metric_pvalues_df[which(cur_metric_pvalues_df$Algorithm != "FCI"), "N"] <- ""

      sink(file=paste0(tables_folder, "pvalues_df_only_valid_", only_valid, "_",
                       cur_metric, "_", data_type, "_sel_top", sel_top, ".txt"))
      print(xtable(cur_metric_pvalues_df,
                   floating=FALSE, latex.environments=NULL,
                   display=display_vec,
                   digits=digits_vec),
            math.style.exponents = TRUE,  include.rownames=FALSE)
      sink()
    }
  }


  ########################
  # Comparing time taken #
  ########################

  time_pvalues_df <- data.frame()
  cur_dcfci_metrics <-  subset(dcfci_metrics, pag_id %in% good_pag_ids)
  all_metric_lists <- list( subset(fci_metrics, pag_id %in% good_pag_ids),
                            subset(cfci_metrics, pag_id %in% good_pag_ids))
  if (data_type == "continuous") {
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(bccd_metrics, pag_id %in% good_pag_ids)
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcd_metrics, pag_id %in% good_pag_ids)
    if ("magsl" %in% methods) {
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(magsl_metrics, pag_id %in% good_pag_ids)
    }
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(gps_metrics, pag_id %in% good_pag_ids)
  }
  all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcfci_metrics, pag_id %in% good_pag_ids)

  metric = "time_taken"
  for (metrics_list_i in 1:(length(methods))) {
    cur_method <- methods[metrics_list_i]
    metrics_list <- all_metric_lists[[metrics_list_i]]
    cur_valid_compat <- 1:nrow(metrics_list)
    #cur_valid_compat <- which(cur_dcfci_metrics$ntop == 1)
    len <- length(cur_valid_compat)
    grpA <- as.numeric(metrics_list[cur_valid_compat, metric])
    if (metrics_list_i <= 3) {
      grpA  <- grpA + 0.04 # time spent with c.i. tests
    }

    grpB <- as.numeric(cur_dcfci_metrics[cur_valid_compat, metric])
    stats_A <- signif(summary(grpA),5)
    stats_B <- signif(summary(grpB),5)
    if (cur_method == "dcfci") {
      pvalue = " "
    } else {
      pvalue <- signif(SignTest(x=grpA, y =grpB, conf.level = 0.95, alternative = "two.sided")$p.value,2)
    }
    decision <- "--"
    if (pvalue < 0.05) {
      if (stats_A[3] < stats_B[3] ||
          stats_A[4] < stats_B[4]) decision = "uparrow"
      else if (stats_A[3] > stats_B[3] ||
               stats_A[4] > stats_B[4]) decision = "downarrow"
    }

    time_pvalues_df <- rbind.data.frame(time_pvalues_df,
                                        data.frame(metric=metric, len=len,
                                                   t(as.numeric(stats_A[1:6])),
                                                   pvalue=pvalue,
                                                   decision
                                        ))
  }

  time_pvalues_df = cbind(
    method=rep(methods_labels[1:(length(methods_labels))], 1), time_pvalues_df)

  time_pvalues_df$method <- factor(time_pvalues_df$method,
                                   levels = methods_labels[1:(length(methods_labels))],
                                   labels = methods_labels[1:(length(methods_labels))])

  time_pvalues_df <- time_pvalues_df[ , -which(colnames(time_pvalues_df) == "metric")]
  colnames(time_pvalues_df) <- c("Algorithm", "nsims",
                                 "Min", "Q1", "Med", "Mean", "Q3", "Max",
                                 "pvalue", "decision")


  sink(file=paste0(tables_folder, "pvalues_df_time_taken", "_",
                   data_type, "_sel_top", sel_top, ".txt"))
  display_vec <- c("s", "s", "d", rep("f", 7), "s")
  digits_vec <- c(0, 0, rep(4, 8), 0)
  print(xtable(time_pvalues_df,
               floating=FALSE, latex.environments=NULL,
               display=display_vec,
               digits=digits_vec),
        math.style.exponents = TRUE,  include.rownames=FALSE)
  sink()


  ###############
  # Checking BIC #
  ################
  #
  # if (data_type == "continuous") {
  #   bic_df <- cbind.data.frame(fci_metrics[, 1:3], fci=fci_metrics$bic,
  #                              cfci=cfci_metrics$bic, bccd=bccd_metrics$bic,
  #                              dcd=dcd_metrics$bic, dcd2=dcd_metrics$dcd_bic,
  #                              magsl=magsl_metrics$bic,
  #                              gps=gps_metrics$bic,
  #                              dcfci=as.numeric(dcfci_metrics$bic))
  #   bic_df <- (bic_df[complete.cases(bic_df),])
  #   dim(bic_df)
  #   agg_bic_df <- aggregate(bic_df[, 4:11],
  #             by = as.list(bic_df[, c("N", "pag_id")]), FUN=mean)
  #   agg_bic_df <- agg_bic_df[order(agg_bic_df$N),]
  #
  #   agg_bic_df
  # }


  ################################
  # Evaluating Difference in BIC #
  # with the True PAG's BIC      #
  # using all PAGs               #
  ################################


  bic_df <- c()
  metric <- "bic"
  for (cur_N in sample_sizes) {

    all_metric_lists <- list( subset(fci_metrics, N == cur_N & pag_id %in% good_pag_ids),
                              subset(cfci_metrics, N == cur_N & pag_id %in% good_pag_ids))
    if (data_type == "continuous") {
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(bccd_metrics, N == cur_N & pag_id %in% good_pag_ids)
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcd_metrics, N == cur_N & pag_id %in% good_pag_ids)
      if ("magsl" %in% methods) {
        all_metric_lists[[length(all_metric_lists) + 1]] <- subset(magsl_metrics, N == cur_N & pag_id %in% good_pag_ids)
      }
      all_metric_lists[[length(all_metric_lists) + 1]] <- subset(gps_metrics, N == cur_N & pag_id %in% good_pag_ids)
    }
    all_metric_lists[[length(all_metric_lists) + 1]] <- subset(dcfci_metrics, N == cur_N & pag_id %in% good_pag_ids)
    cur_true_bic <- as.numeric(subset(true_bic, N == cur_N & pag_id %in% good_pag_ids)[, "bic", drop=TRUE])

    for (metrics_list_i in 1:length(all_metric_lists)) {
      metrics_list <- all_metric_lists[[metrics_list_i]]
      cur_method <- methods[[metrics_list_i]]
      method_bic <- cbind(method=cur_method, metrics_list[, c("N", "pag_id", metric)])
      method_bic <- cbind(method_bic, true=cur_true_bic)
      bic_df <- rbind(bic_df, method_bic)
    }
  }

  bic_df$method <- factor(bic_df$method,
                             levels = methods,
                             labels = methods_labels,
                          ordered = TRUE)
  bic_df$N <- factor(bic_df$N, ordered = TRUE)
  bic_df$pag_id <- as.factor(bic_df$pag_id)
  bic_df$diff_bic <- bic_df$bic - bic_df$true

  min_bic <- floor(min(bic_df$diff_bic, na.rm=T))
  max_bic <- ceiling(max(bic_df$diff_bic, na.rm=T))


  bic_df$diff_bic2 <- cut(bic_df$diff_bic,
                              breaks=c(min_bic, -5000,
                                       #-1000,
                                       -100,
                                       -0.0001, 0.0001,
                                       100,
                                       #1000,
                                       5000, max_bic),
                          labels=c("< -5000",
                                   #"[-5000, 0)",
                                   "[-5000, -100)", "[-100,0)",
                                   "0",
                                   "(0, 100]", "(100, 5000]",
                                   #"(0, 5000]",
                                   "> 5000"),
                          ordered_result = TRUE)
  table(bic_df$diff_bic2)

  table(subset(bic_df, method == "dcFCI" & N == 50000)$diff_bic2)
  table(subset(bic_df, method == "FCI" & N == 50000)$diff_bic2)
  table(subset(bic_df, method == "cFCI" & N == 50000)$diff_bic2)


  agg_bic_diff_df <- data.frame()
  for (cur_method in methods_labels) {
    for (cur_N in sample_sizes) {
      cur_bic_diff_df <- subset(bic_df, N == cur_N & method == cur_method)
      cur_agg_bic_diff_df <- cbind.data.frame(diff_bic = levels(cur_bic_diff_df$diff_bic2),
                                              count = as.numeric(table(cur_bic_diff_df$diff_bic2)),
                                              Algorithm = cur_method,
                                              N = cur_N)
      agg_bic_diff_df <- rbind.data.frame(agg_bic_diff_df, cur_agg_bic_diff_df)
    }
  }


  agg_bic_diff_df$diff_bic <- factor(agg_bic_diff_df$diff_bic,
                                     levels=levels(bic_df$diff_bic2),
                                     labels=levels(bic_df$diff_bic2))
  agg_bic_diff_df$Algorithm <- factor(agg_bic_diff_df$Algorithm,
                                      levels = levels(bic_df$method),
                                      labels=levels(bic_df$method))

  agg_bic_diff_df$N <- factor(agg_bic_diff_df$N, ordered = TRUE)

  labels_bic_diff <- levels(agg_bic_diff_df$diff_bic)

  ggplot(agg_bic_diff_df, aes(x=diff_bic, y=count, fill = Algorithm)) +
    geom_col(position = position_dodge(width=0.8), width=0.8) +
    #geom_col(position = "dodge", width=0.5) +
    #geom_bar(stat="identity") +
    ylab("Count") +
    ggtitle("Difference in BIC by Algorithm and Sample Size") +
    scale_x_discrete("Difference in BIC (Estimated - True)", labels=labels_bic_diff) +
    scale_fill_manual(values = method_colors) +
    facet_wrap(~ N, ncol=1) +
    theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5, size=14))

  ggsave(paste0(figures_folder, "diff_bin", "_", data_type, "_sel_top_", sel_top, ".png"),
         width=8, height =6)
}
}


