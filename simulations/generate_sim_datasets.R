rm(list=ls())

library(FCI.Utils)
library(dcFCI)
library(jsonlite)

library(doFuture)
library(future.apply)
library(rje)
library(dplyr)
source("./simulations/sim_utils.R")

n_cores <- 12
plan("multicore", workers = n_cores) # forking
# plan("cluster", workers = n_cores)


sim_true_pags_file <- paste0("./simulations/random_p5_true_pags_list.RData")
load(sim_true_pags_file)
#lapply(true_pags_list, renderAG)


sample_sizes = c(200, 500, 1000, 5000, 10000, 50000)
pag_ids <- 1:length(true_pags_list)
sim_ids = 1:30

#######################
# Generating datasets #
#######################

runMAGSL = TRUE

for (data_type in c("continuous", "mixed")) {
  for (pag_id in pag_ids) {
    for (N in sample_sizes) {
      for (sim in sim_ids) {

        cat("PAG: ", pag_id, "N: ", N, "sim: ", sim, "\n")


        dat_folder <- paste0("../dcFCI_Simulations/",
                             data_type, "/",  sprintf("%02d", pag_id),
                             "/n", N, "/sim_", sim, "/")

        if (!file.exists(dat_folder)) {
          dir.create(dat_folder, recursive = TRUE)
        }

        true.amat.pag <- true_pags_list[[pag_id]]

        dat_file <- list.files(dat_folder, pattern=glob2rx("dat_*.csv"))

        if (length(dat_file) == 0) {
          f_dat_out <- generateDatasetHelper(true.amat.pag, data_type, N,
                                            max_tries = 1000)
          dat <- f_dat_out$dat
          suffStat <- f_dat_out$suffStat
          cur_seed <- f_dat_out$seed

          dat_file <- paste0("dat_seed_", cur_seed, ".csv")
          if (!is.null(dat)) {
            write.csv(dat, file=paste0(dat_folder, dat_file), row.names = FALSE)

            if (runMAGSL) {
              magsl_csv <- paste0(dat_folder, "/", "magsl_", dat_file)
              write.table(dat, file=magsl_csv,
                            row.names = FALSE, col.names = FALSE, sep = " ")
            }

            save(suffStat, file=paste0(dat_folder, "suffStat_sim_", sim, ".RData"))
          } else {
            message("Problem while generating dataset.\n")
            next()
          }
        }
      }
    }
  }
}
