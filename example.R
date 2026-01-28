rm(list=ls())

library(FCI.Utils)
library(dcFCI)

library(jsonlite)
library(rje)
library(dplyr)

# If running in parallel:
run_parallel = TRUE

if (run_parallel) {
  require(doFuture)
  require(future.apply)
  n_cores <- 2 # change to the number of cores/CPU available
  #plan("multisession", workers = n_cores) # to run on RStudio
  plan("multicore", workers = n_cores) # using forking on RStudio
  #plan("cluster", workers = n_cores) # to run on a server / terminal
}

#####################
#  Paper's Figure 1 #
#####################

# Generating true ADMG and PAG
g_type <- "discr1_nc"
true.admg <- getDAG(g_type)
true.amat.pag <- getTruePAG(true.admg$dagg, verbose = FALSE)@amat
true.sepset <- getPAGImpliedSepset(true.amat.pag)

renderAG(true.amat.pag, add_index = F)
formatSepset(true.sepset)

## Generating dataset ##
########################


run_paper_example = TRUE
if (run_paper_example) {
  N <- 10000
  data_type = "continuous"
  seed = 1042602069
  set.seed(seed)
  b.default <- sapply(1:nrow(dagitty::edges(true.admg$dagg)),
                    function(x) {runif(1, -0.6, 0.6)})
  dat <- generateDataset(adag = true.admg$dagg, N=N,
                             type=data_type, b.default = b.default)$dat

  # Alternatively:
  # dat_folder <-"./data_example/"
  # dat_file <- "dat_seed_1042602069.csv"
  # dat <- read.csv(file=paste0(dat_folder, dat_file))
} else {
  # Generating a random dataset according to the true ADMG

  # data_type = "continuous" # for only continuous (Gaussian) variables
  # data_type = "mixed" # for a mixture of continuous and discrete variables
  data_type = "continuous"
  N <- 1000 # sample size

  # For interesting examples with N == 1000:
  selected_seed = TRUE
  if (selected_seed && N == 1000 && data_type == "mixed") {
    seed <- 735706305
  } else if (selected_seed && N == 1000 && data_type == "continuous") {
    seed <- 47555208
  } else {
    seed <- sample(1:.Machine$integer.max, 1)
  }

  dat <- NULL
  while (is.null(dat)) {
    if (!selected_seed) {
      seed <- sample(1:.Machine$integer.max, 1)
    }
    set.seed(seed)

    f.args <- NULL
    if (data_type == "mixed") {
      all_vars = names(dagitty::coordinates(true.admg$dagg)[["x"]])
      f.args <- list()

      # 1: continuous; 2: binary; 3: multinominal
      # var_levels <- sample(1:3, size=length(all_vars), replace = TRUE)
      var_levels <- c(1,2,1,2,3) # for "A", "B", "Uab", "X", "Y"

      for (vari in 1:length(all_vars)) {
        var_name <- all_vars[vari]
        f.args[[var_name]] <- list(levels = var_levels[vari])
      }
    }
    dat_out <- generateDataset(adag = true.admg$dagg, N=N,
                               type=data_type, f.args=f.args)
    dat <- dat_out$dat
  }
}


## Computing conditional independence tests in parallel ##
##########################################################

vars_names <- colnames(dat)
covs_names = c()
indepTest <- mixedCITest

suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
vars_df <- dat[,vars_names, drop=FALSE]

citestResults <- getAllCITestResults(
  vars_df, indepTest, suffStat, m.max = Inf, # limit m.max for larger graphs
  computeProbs = TRUE,        # Required if running dcFCI
  eff_size = 0.05,            # Minimum effect size threshold when computeProbs = TRUE
  fileid = paste("example_seed_", seed),  # Unique identifier for the dataset / analysis
  recover_citestResults = TRUE, # Set TRUE to resume computation after a crash
  results_folder= file.path(getwd(), "tmp","citestsResults"), # same as default -- change if needed
  log_folder = file.path(getwd(), "tmp", "logs")) # same as default -- change if needed


suffStat$citestResults <- citestResults # so we don´t redo tests later

# Checking true PAG's scores

# Straightforward PAG score: comparable, as it includes all tests
true_sf_score <- getStraightforwardPAGScore2(true.amat.pag, suffStat)
true_sf_score

# MEC-targeted PAG score
true_mec_score <- getMECTargetedPAGScore(true.amat.pag, suffStat, ord=Inf)
true_mec_score


##########################################
# Setting up causal discovery parameters #
##########################################

m.max = Inf
labels <- vars_names
p <- length(labels)
fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
fixedGaps = NULL
alpha = 0.01
verbose = TRUE

###############
# Running FCI #
###############

# run original FCI
conservative = FALSE # change to TRUE for conservative FCI (cFCI)
maj.rule = FALSE

start_time <- Sys.time()
fit_fci <- pcalg::fci(suffStat, indepTest = mixedCITest,
                      skel.method = "stable", labels = labels, m.max=m.max,
                      NAdelete = NAdelete, type = "normal", alpha = alpha,
                      verbose = verbose, conservative = conservative,
                      maj.rule = maj.rule)
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

fci_pag <- fit_fci@amat
# renderAG(fci_pag, add_index = FALSE)

fci_metrics <- getMetrics(true.amat.pag, fci_pag,
                          est.sepset=getPAGImpliedSepset(fci_pag),
                          dat=NULL,  # specify dat if BIC should be computed
                          conservative = FALSE)

fci_sepset <- fixSepsetList(fit_fci@sepset)
formatSepset(fci_sepset)

fci_metrics

#################
# Running dcFCI #
#################

sel_top = 1 # parameter k in the paper


# alpha=0.05; m.max = Inf; fixedGaps = NULL; fixedEdges = NULL;
# verbose = TRUE; sel_top = 1; prob_sel_top = FALSE;
# run_parallel = TRUE; allowNewTests=FALSE;
# list.max = 500; pH0Thresh=0.8;
# log_folder = file.path(getwd(), "tmp", "logs")

start_time <- Sys.time()
fit_dcfci <- dcFCI(suffStat, indepTest, labels, alpha,
                   m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                   verbose = FALSE, # 1+: how verbose it can be
                   sel_top = 1,
                   run_parallel = TRUE,
                   allowNewTests=TRUE,
                   list.max = 1500, # may be increased for highly uncertain cases
                   pH0ThreshMin=0.3, # if pH0ThreshMin <= pH0 <= pH0ThreshMax, it will be a potential conditional independence
                   pH0ThreshMax = 1,
                   log_folder = file.path(getwd(), "tmp", "logs"))
end_time <- Sys.time()

time_taken <- end_time - start_time
cat("\n\n dcFCI runtime: ", time_taken, "\n\n")


length(fit_dcfci$top_dcPAGs)

top_dcfci_pag <- fit_dcfci$allPAGList[[1]]$amat.pag
renderAG(top_dcfci_pag, add_index = FALSE)

top_dcfci_sepset <- fit_dcfci$allPAGList[[1]]$sepset
formatSepset(top_dcfci_sepset)

dcfci_metrics <- getMetrics(true.amat.pag, top_dcfci_pag,
                            est.sepset=top_dcfci_sepset,
                            dat=NULL,  # specify dat if BIC should be computed
                            conservative = FALSE)
dcfci_metrics


# MEC-targeted PAG scores:
validPAGScores <- subset(fit_dcfci$mec_score_df,
                         violations == FALSE & duplicated == FALSE)

head(validPAGScores, 10) # Top 10 PAGs
#validPAGScores # All evaluated PAGs

# sepsets that led to valid PAGs, ranked by MEC-Targeted PAG score
#lapply(fit_dcfci$allPAGList[validPAGScores$pag_list_id],
#       function(x) {formatSepset(x$sepset)})



########################
# Printing the results #
########################

renderAG(fci_pag, add_index = F) # FCI PAG
cat("FCI PAG's SHD: ", fci_metrics$shd, "\n")

renderAG(fit_dcfci$top_dcPAGs[[1]]$amat.pag, add_index = F) # dcFCI PAG
cat("dcFCI PAG's SHD: ", dcfci_metrics$shd, "\n")



