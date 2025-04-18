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
true.amat.pag <- getTruePAG(true.admg$dagg, verbose = FALSE)@amat
renderAG(true.amat.pag, add_index = F)


# Generating a dataset according to the true ADMG
data_type = "mixed" #continuous
N <- 10000

if (data_type == "continuous") {
  set.seed(1042602069)
  b.default <- sapply(1:nrow(dagitty::edges(true.admg$dagg)),
                      function(x) {runif(1, -0.6, 0.6)})
  dat_out <- generateDataset(adag = true.admg$dagg, N=N,
                             type=data_type, b.default = b.default)
} else if (data_type == "mixed") {
  all_vars = names(dagitty::coordinates(true.admg$dagg)[["x"]])
  f.args <- list()

  #cur_seed <- sample(1:.Machine$integer.max, 1)
  #set.seed(cur_seed)
  set.seed(1412620727)

  # 1: continuous; 2: binary; 3: multinominal
  # var_levels <- sample(1:3, size=length(all_vars), replace = TRUE)

  var_levels <- c(1,2,1,2,3) # for "A", "B", "Uab", "X", "Y"

  for (vari in 1:length(all_vars)) {
    var_name <- all_vars[vari]
    f.args[[var_name]] <- list(levels = var_levels[vari])
  }

  dat_out <- generateDataset(adag = true.admg$dagg, N=N,
                                     type=data_type, f.args=f.args,
                                     coef_thresh = 0.1)
}

dat <- dat_out$dat

# Computing conditional independence tests in parallel

vars_names <- colnames(dat)
covs_names = c()
indepTest <- mixedCITest

suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
vars_df <- dat[,vars_names, drop=FALSE]
citestResults <- getAllCITestResults(vars_df, indepTest, suffStat,
                                     m.max=Inf, computeProbs = TRUE,
                                     eff_size = 0.01 # ensure a mininum effect size
                                     #eff_size = NULL # uses max RMSE
                                     )

suffStat$citestResults <- citestResults # so we don´t redo tests later

# Checking faithfulness degree
faithf_degree <- getFaithfulnessDegree(true.amat.pag, citestResults, alpha = 0.01)
faithf_degree$f_citestResults


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

fci_pag <- fit_fci@amat

fci_metrics <- getMetrics(true.amat.pag, fci_pag,
                          est.sepset=getPAGImpliedSepset(fci_pag),
                          dat=NULL,  # specify dat if BIC should be computed
                          conservative = FALSE)

fci_sepset <- fixSepsetList(fit_fci@sepset)
formatSepset(fci_sepset)

fci_metrics$shd

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

dcfci_pag <- fit_dcfci$top_dcPAGs[[1]]$amat.pag


dcfci_metrics <- getMetrics(true.amat.pag, dcfci_pag,
                          est.sepset=getPAGImpliedSepset(dcfci_pag),
                          dat=NULL,  # specify dat if BIC should be computed
                          conservative = FALSE)

dcfci_sepset <- fit_dcfci$top_dcPAGs[[1]]$sepset
formatSepset(dcfci_sepset)
fci_metrics$shd


########################
# Printing the results #
########################

renderAG(fci_pag, add_index = F) # FCI PAG
cat("FCI PAG's SHD: ", fci_metrics$shd, "\n")

renderAG(fit_dcfci$top_dcPAGs[[1]]$amat.pag, add_index = F) # dcFCI PAG
cat("dcFCI PAG's SHD: ", dcfci_metrics$shd, "\n")



