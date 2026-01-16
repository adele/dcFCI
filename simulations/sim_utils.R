generateDatasetHelper <- function(true.amat.pag, data_type, N,
                                  score_upp_thresh = 0.05,
                                  eff_size = 0.05,
                                  m.max=Inf, max_tries=20,
                                  verbose=FALSE) {

  stats <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(stats) <- c("seed", "frechetUB", "mse")

  old_var_names <- colnames(true.amat.pag)

  if (data_type == "mixed") {
     score_upp_thresh = 0.005
  }

  ntries = 0
  while (ntries < max_tries) {
    ntries = ntries + 1

    if (verbose)
      message(paste0("Generating sample -- ntries: ", ntries))

    cur_seed <- sample(1:.Machine$integer.max, 1)
    set.seed(cur_seed)

    colnames(true.amat.pag) <- paste0("V", 1:ncol(true.amat.pag)) #c("W", "X", "Y", "Z") #
    rownames(true.amat.pag) <- colnames(true.amat.pag)
    new_labels <- colnames(true.amat.pag)

    f.args <- NULL
    if (data_type == "mixed") {
      f.args <- list()
      if (length(new_labels) == 4) {
        var_levels <- c(1, 3, 2, 1)
      } else if (length(new_labels) == 5) {
        var_levels <- c(1, 1, 1, 3, 2)
      }

      for (vari in 1:length(new_labels)) {
        var_name <- new_labels[vari]
        f.args[[var_name]] <- list(levels = var_levels[vari])
      }
    }

    # Generating the dataset with variables as columns and observations as rows
    #adat_out <- generateDataset(adag = adag_out$dagg, N=N, type=data_type)
    adat_out <- generateDatasetFromPAG(apag = true.amat.pag, N=N,
                                       type=data_type, f.args=f.args,
                                       coef_thresh = 0.2, verbose=TRUE)
    dat <- adat_out$dat
    if (is.null(dat)) {
      if (verbose)
        message("Cannot created dat for such a PAG.")
      next
    }

    if (verbose)
      message("Dataset created for such a PAG.")

    colnames(true.amat.pag) <- old_var_names
    rownames(true.amat.pag) <- old_var_names
    colnames(dat) <- old_var_names

    cat_dat <- dat[, sapply(dat, is.factor)]
    if (!all(sapply(cat_dat, function(x) { all(table(x) > 0) } ))) {
      next
    }
    #str(dat)

    ############################################
    # Computing Conditional Independence Tests #
    ############################################

    vars_names <- colnames(true.amat.pag) # only observed variables
    covs_names = c()
    indepTest <- mixedCITest


    suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
    vars_df <- dat[,vars_names, drop=FALSE]
    citestResults <- getAllCITestResults(vars_df, indepTest, suffStat,
                                         m.max=m.max+1, computeProbs = TRUE,
                                         eff_size = eff_size)
    suffStat$citestResults <- citestResults

    # sf_score2 <- getStraightforwardPAGScore(true.amat.pag, citestResults, m.max)
    # sf_score2$min_dep_test
    # sf_score2$min_indep_test

    sf_score <- getStraightforwardPAGScore2(true.amat.pag, suffStat, m.max)
    stats <- rbind.data.frame(stats,
      data.frame(seed=cur_seed, frechetUB=sf_score$frechetUB, mse=sf_score$mse))

    if (sf_score$frechetUB > score_upp_thresh) {
      return(list(seed=cur_seed, dat=dat, suffStat = suffStat, stats=stats))
    } else {
      if (verbose)
        print("Too unfaithful dataset... ")
    }
  }
  return(list(seed=NULL, dat=NULL, suffStat = NULL, stats=stats))
}

# TODO: Remove after compiling FCI.Utils

isValidAG <- function(amat.ag) {
  isAGret <- tryCatch({
    if (!is.null(amat.ag)) {
      ug_ag <- (amat.ag == 3 & t(amat.ag == 3)) * 1
      bg_ag <- (amat.ag == 2 & t(amat.ag == 2)) * 1
      dg_ag <- (amat.ag == 2 & t(amat.ag == 3)) * 1
      ag_ggm <- ggm::makeMG(dg_ag, ug_ag, bg_ag)
      ggm::isAG(ag_ggm)
    } else {
      FALSE
    }
  },
  error=function(cond) {
    print(cond)
    return(FALSE)
  },
  warning=function(cond) {
    print(cond)
    return(FALSE)
  })
}

