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
  n_cores <- 12
  #plan("multisession", workers = n_cores)
  plan("multicore", workers = n_cores) # forking
}

#####################
#  Paper's Figure 1 #
#####################

# Generating true ADMG and PAG
g_type <- "discr1_nc"
true.admg <- getDAG(g_type)
true.amat.pag <- getTruePAG(true.admg$dagg, verbose = TRUE)@amat
renderAG(true.amat.pag, add_index = F)


# Generating a dataset according to the true ADMG
data_type = "continuous"
N <- 10000

set.seed(1042602069)
dat_out <- generateDataset(adag = true.admg$dagg, N=N, type=data_type)
dat <- dat_out$dat


# Computing conditional independence tests in parallel

vars_names <- colnames(dat)
covs_names = c()
indepTest <- mixedCITest

suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
vars_df <- dat[,vars_names, drop=FALSE]
citestResults <- getAllCITestResults(vars_df, indepTest, suffStat,
                                     m.max=Inf, computeProbs = TRUE)
suffStat$citestResults <- citestResults # so we don´t redo tests later

# Checking faithfulness degree
faithf_degree <- getFaithfulnessDegree(true.amat.pag, citestResults, alpha = 0.01)



##########################################
# Setting up causal discovery parameters #
##########################################

m.max = Inf
labels <- vars_names
p <- length(labels)
fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
fixedGaps = NULL
NAdelete = FALSE
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


renderAG(fit_fci@amat, add_index = F)




#################
# Running dcFCI #
#################

sel_top = 1 # parameter k in the paper

start_time <- Sys.time()
fit_dcfci <- dcFCI(suffStat, mixedCITest,
                   labels=labels, alpha=alpha,
                   NAdelete = NAdelete,  m.max=m.max,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   verbose = verbose, sel_top = sel_top,
                   run_parallel = run_parallel)
end_time <- Sys.time()
time_taken <- end_time - start_time

# MEC-targeted PAG scores in each iteration:
fit_dcfci$top_scoresDF

# MEC-targeted PAG scores in each iteration:
ntop <- length(fit_dcfci$top_dcPAGs)
ntop # number of PAGs ranked at 1st.

renderAG(fit_dcfci$top_dcPAGs[[1]]$amat.pag, add_index = F)


