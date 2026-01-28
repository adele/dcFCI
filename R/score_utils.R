# library(lavaan)

#' @importFrom lavaan summary
#' @importFrom FCI.Utils getMAG getIgraphMAG dagitty2amat
#' @importFrom SEMgraph SEMrun
#' @importFrom dagitty toMAG
#' @importFrom stats shapiro.test
#' @importFrom RNOmni RankNorm
#' @importFrom methods is
#' @export getBIC
getBIC <- function(amat.pag, dat, gaussian_vars=NULL, group_var=NULL,
                   type="adjBIC", alpha=0.05, applyINT=FALSE) {
  if (is.null(gaussian_vars)) {
    gaussian_vars <- colnames(amat.pag)
  }

  cur_mag <- getMAG(amat.pag)
  can_dag <- dagitty::canonicalize(cur_mag$magg)$g
  hidden_vars <- setdiff(colnames(amat.pag), c(gaussian_vars, group_var))
  dagitty::latents(can_dag) <- c(dagitty::latents(can_dag), hidden_vars)
  proj_mag <- dagitty::toMAG(can_dag)
  ig.mag <- getIgraphMAG(dagitty2amat(proj_mag))
  #plot(ig.mag)
  #sel_dat <- adat[,dagitty::names.dagitty(proj_mag)]

  sel_dat <- dat[,c(group_var, gaussian_vars)]
  if (!is.null(group_var) && methods::is(sel_dat[,group_var], "factor")) {
    sel_dat[,group_var] <- as.numeric(as.character(sel_dat[,group_var]))
  }

  for(gvar in gaussian_vars) {
    if (applyINT) {
      shp <- stats::shapiro.test(sel_dat[, gvar])
      if (shp$p.value < alpha) {
        cat("WARNNIG: variable ", gvar,
            "does not follow a Gaussian distribution; p.value=", shp$p.value, "\n")
        cat("    --> applying inverse normal transformation.. \n")
        sel_dat[, gvar] <- RNOmni::RankNorm(sel_dat[, gvar])
        shp <- stats::shapiro.test(sel_dat[, gvar])
        if (shp$p.value < alpha) {
          cat("WARNNIG: variable ", gvar,
              "still does not follow a Gaussian distribution; p.value=", shp$p.value, "\n")
        } else {
          cat("Variable", gvar, "now follows a Gaussian distribution; p.value=",
              shp$p.value, "\n")
        }
      }
    }
  }

  fit <- tryCatch({
    if (!is.null(group_var)) {
      SEMgraph::SEMrun(graph=ig.mag, data=sel_dat, group=sel_dat[,group_var]) #, algo = "ricf")
    } else {
      SEMgraph::SEMrun(graph=ig.mag, data=sel_dat, group=NULL) #, algo = "ricf")
    }
  }, error = function(cond) {
    message("Error computing BIC")
    message(conditionMessage(cond))
    NA
  })

  if (length(fit) == 1 && is.na(fit)) {
    return(list(bic=NA, fit=NA, ig.mag=ig.mag))
  } else {
    summ_fit <- lavaan::summary(fit$fit, fit.measures = TRUE)
    if (type == "adjBIC") {
      return(list(bic=summ_fit$fit["bic2"], fit=fit, ig.mag=ig.mag))
    } else {
      return(list(bic=summ_fit$fit["bic"], fit=fit, ig.mag=ig.mag))
    }
  }
}

#' @export getStraightforwardPAGScore
getStraightforwardPAGScore <- function(apag, citestResults, m.max = Inf) {
  if (is.null(citestResults)) {
    cat("Error: citestResults must contain all required CI test results.\n")
    return(NULL)
  }

  na_ids <- which(is.na(citestResults$pvalue))
  if (length(na_ids) > 0) {
    if (length(na_ids) > 0) {
      citestResults$pvalue[na_ids] <- NA
      citestResults$pH0[na_ids] <- 0.5
      citestResults$pH1[na_ids] <- 0.5
    }
  }

  faithf <- getFaithfulnessDegree(apag, citestResults)
  indep_tests <- subset(faithf$f_citestResults, type == "indep")
  if (!is.infinite(m.max)) {
    indep_tests <- subset(faithf$f_citestResults, type == "indep" & ord <= m.max+1)
  }
  indep_probs <- indep_tests$pH0
  min_indep_test <- indep_tests[which.min(indep_probs),]

  dep_tests <- subset(faithf$f_citestResults, type == "dep")
  if (!is.infinite(m.max)) {
    dep_tests <- subset(faithf$f_citestResults, type == "dep" & ord <= m.max+1)
  }
  dep_probs <- dep_tests$pH1
  min_dep_test <- dep_tests[which.min(dep_probs),]

  score <- dcFCI::getProbConjunction(c(indep_probs, dep_probs))
  mse <- dcFCI::getNormalizedSquaredL2DistanceFromCertainty(c(indep_probs, dep_probs))

  return(list(score=score,
              mse=mse,
              min_indep_test = min_indep_test,
              min_dep_test = min_dep_test,
              indep_tests = indep_tests,
              false_indep_tests = subset(indep_tests, bf==FALSE),
              dep_tests = dep_tests,
              false_dep_tests = subset(dep_tests, bf==FALSE)))
}

#' @export getProbConjunctionFromIntervals
getProbConjunctionFromIntervals <- function(prob_intervals) {
  min_int <- getProbConjunction(sapply(prob_intervals, function(x) {x[1]}))
  max_int <- getProbConjunction(sapply(prob_intervals, function(x) {x[2]}))
  vals <- c(min_int, max_int)
  return(c(min(vals), max(vals)))
}

# equivalences:
# sum(probs) - (length(probs)-1) --> Frechet's lower bound
# 1 - sum(1 - probs) --> 1 - (maximum prob of at least one event fails).
# Note that (1 - prob_i) is the is the maximum possible probability that at
# least one event fails, given only the marginals. Equality holds only when the
# failure events are disjoint.
# the probability that all events succeed is  1 - sum(1 - probs), i.e, at least
# one minus the maximum possible probability that some event fails
# sum(1 - probs) is the total slack away from certainty across all events.
# sum(1 - probs) is also the l1 distance from 1.
#' @export getL2DistanceFromCertainty
getL2DistanceFromCertainty <- function(probs) {
  l2d <- sqrt(sum((rep(1, length(probs)) - probs)^2))
  return(l2d)
}

#' @export getSquaredL2DistanceFromCertainty
getSquaredL2DistanceFromCertainty <- function(probs) {
  l2d2 <- sum((rep(1, length(probs)) - probs)^2)
  return(l2d2)
}

# Same as MSE
#' @export getNormalizedSquaredL2DistanceFromCertainty
getNormalizedSquaredL2DistanceFromCertainty <- function(probs) {
  mse <- getSquaredL2DistanceFromCertainty(probs)/length(probs)
  return(mse)
}

#' @export getProbConjunction
getProbConjunction <- function(probs) {
  minProb <- max(0, sum(probs) - (length(probs)-1))
  maxProb <- min(probs)
  return(c(minProb, maxProb))
}

#' @export getProbDisjuntionFromIntervals
getProbDisjuntionFromIntervals <- function(prob_intervals) {
  min_int <- getProbDisjunction(sapply(prob_intervals, function(x) {x[1]}))
  max_int <- getProbDisjunction(sapply(prob_intervals, function(x) {x[2]}))
  vals <- c(min_int, max_int)
  return(c(min(vals), max(vals)))
}

#' @export getProbDisjunction
getProbDisjunction <- function(probs) {
  minProb <- max(probs)
  maxProb <- min(1, sum(probs))
  return(c(minProb, maxProb))
}


# pH0 for one or more conditional independence
# tests corresponding to implied independencies.
# It returns an score for the probability of all (conjunction)
# independencies being true.
#' @export getIndepScore
getIndepScore <- function(h0_probs, ret_int=TRUE) {
  indep_score_int <- getProbConjunction(probs=h0_probs)
  if (ret_int) {
    return(indep_score_int)
  } else {
    return(indep_score_int[1]) # the lower bound
  }
}

#' @export getDepScore
getDepScore <- function(h1_probs, ret_int=TRUE) {
  indep_score_int <- getProbConjunction(probs=h1_probs)
  if (ret_int) {
    return(indep_score_int)
  } else {
    return(indep_score_int[1]) # the lower bound
  }
}

#' @export getCITestProbs
getCITestProbs <- function(citestResults, n) {
  #TODO fix warning
  probs <- foreach (i = 1:nrow(citestResults), #.options.snow = opts,
             .verbose = TRUE, .export = ls(globalenv())) %dopar% {
               return(unlist(pvalue2probs(citestResults[i, "pvalue"], n=n)))
             }
  probs <- sapply(probs, function(x) {x})
  citestResults <- cbind.data.frame(citestResults, t(probs))
  return(citestResults)
}

getIntervalProduct <- function(a, b) {
  prods <- c(a[1] * b[1], a[1] * b[2], a[2] * b[1], a[2] * b[2])
  return(c(min(prods), max(prods)))
}

getTrivialScore <- function(amat.pag, citestResults, ord, verbose=FALSE) {
    labels <- colnames(amat.pag)

    ord_citestResults <- subset(citestResults, ord <= ord)
    f_citestResults <- c()
    for (i in 1:nrow(ord_citestResults)) {
      cur_row <- ord_citestResults[i, , drop=TRUE]
      snames <- labels[FCI.Utils::getSepVector(cur_row$S)]
      xname <- labels[cur_row$X]
      yname <- labels[cur_row$Y]

      def_msep <- FCI.Utils::isMSeparated(amat.pag, xname, yname, snames,
                               verbose=verbose)
      if (def_msep) {
        ret <- c(cur_row, type="indep")
        f_citestResults <- rbind.data.frame(f_citestResults, ret)
      } else {
        ret <- c(cur_row, type="dep")
        f_citestResults <- rbind.data.frame(f_citestResults, ret)
      }
    }

    score <- getProbConjunction(c(subset(f_citestResults, type == "indep")$pH0,
      subset(f_citestResults, type == "dep")$pH1))
    return(score)
}

getAnnotatedCITestResults <- function(amat.pag, ord=NULL, citestResults, verbose=FALSE) {
  labels <- colnames(amat.pag)

  if (is.null(ord)) {
    ord <- ncol(amat.pag) - 2
  }

  a_citestResults <- c()
  for (i in 1:nrow(citestResults)) {
    cur_row <- citestResults[i, , drop=TRUE]
    snames <- labels[FCI.Utils::getSepVector(cur_row$S)]
    xname <- labels[cur_row$X]
    yname <- labels[cur_row$Y]

    if (amat.pag[cur_row$X, cur_row$Y] != 0 &&
        amat.pag[cur_row$Y, cur_row$X] != 0 &&
        cur_row$ord > ord) {
      # there is an edge and ord has not been evaluated yet
      ret <- c(cur_row, type="dep")
      if (ret$pH1 < 0.5) {
        ret$pH0 <- ret$pH1 <- 0.5 # giving unknown probability to it
      }
      a_citestResults <- rbind.data.frame(a_citestResults, ret)
    } else {
      def_msep <- FCI.Utils::isMSeparated(amat.pag, xname, yname, snames,
                               verbose=verbose)
      if (def_msep) {
        ret <- c(cur_row, type="indep")
        a_citestResults <- rbind.data.frame(a_citestResults, ret)
      } else {
        ret <- c(cur_row, type="dep")
        a_citestResults <- rbind.data.frame(a_citestResults, ret)
      }
    }
  }
  return(a_citestResults)
}


# This function computes a score based on the different citestResults
# of each candidate PAG in cur_ord_pag_list, considering the current list ``ord''
# for invalid candidate PAGs, score is [0,0]
# in case of no difference exists among the candidates, score is [1,1]
getSymmDiffCITestResults <- function(cur_ord_pag_list, ord,
                                     verbose=FALSE) {
  valid_pags_ids <- which(!sapply(cur_ord_pag_list, function(x) { x$violations }))

  if (length(valid_pags_ids) == 0) {
    return(NULL)
  }
  cur_ord_valid_pags <- cur_ord_pag_list[valid_pags_ids]
  eval_citests_List <- lapply(cur_ord_valid_pags, function(x) { x$mec$all_citests} )
  all_eval_citests <- bind_rows(eval_citests_List)
  unique_ci_relations <- all_eval_citests[!duplicated(all_eval_citests[,1:4]), ]

  new_eval_citests_List <- list()
  for (j in seq_along(1:length(eval_citests_List))) {
    new_eval_citests_List[[length(new_eval_citests_List) + 1]] <- data.frame()
    if (nrow(unique_ci_relations) > 0) {
      for (i in 1:nrow(unique_ci_relations)) {
        cur_ci_relation = unique_ci_relations[i, 1:7]
        cur_ci_info <- subset(eval_citests_List[[j]], X == cur_ci_relation$X &
                              Y == cur_ci_relation$Y &
                              S == cur_ci_relation$S)
        if (nrow(cur_ci_info) == 0) {
          labels <- colnames(cur_ord_valid_pags[[j]]$amat.pag)

          snames <- labels[FCI.Utils::getSepVector(cur_ci_relation$S)]
          xname <- labels[cur_ci_relation$X]
          yname <- labels[cur_ci_relation$Y]
          def_msep <- FCI.Utils::isMSeparated(
            cur_ord_valid_pags[[j]]$amat.pag, xname, yname, snames,
                                              verbose=verbose)
          if (def_msep) {
            ret <- c(cur_ci_relation, type="indep", structure="new")
            new_eval_citests_List[[j]] <- rbind.data.frame(
              new_eval_citests_List[[j]], ret)
          } else {
            ret <- c(cur_ci_relation, type="dep", structure="new")
            new_eval_citests_List[[j]] <- rbind.data.frame(
              new_eval_citests_List[[j]], ret)
          }
        } else {
          new_eval_citests_List[[j]] <- rbind.data.frame(
            new_eval_citests_List[[j]], cur_ci_info)
        }
      }
    }
  }

  symmDiffCITests <- data.frame()
  symm_diff_ids <- c()
  if (nrow(unique_ci_relations) > 0 && length(new_eval_citests_List) > 1) {
    deps_indeps <- sapply(new_eval_citests_List, function(x) { x$type} )

    symm_diff_ids <- which(apply(deps_indeps, 1, function(x) {
        length(unique(x)) > 1
      }))

    symmDiffCITests <- unique_ci_relations[symm_diff_ids, 1:7]
  }

  diff_citestResults_list <- list()
  diff_citestResults_list[1:length(cur_ord_pag_list)] <- list(data.frame())
  diff_citestResults_list[valid_pags_ids] <-
    lapply(new_eval_citests_List, function(x) {
      if (length(symm_diff_ids) > 0) {
        x[symm_diff_ids, ]
      } else {
        data.frame()
      }
    })

  #lapply(diff_citestResults_list[valid_pags_ids], function(x) { dim(x) })


  for (i in seq_len(length(cur_ord_pag_list))) {
    cur_diff_citestResults <- diff_citestResults_list[[i]]
    cur_ord_pag_list[[i]]$diff_ord_citestResults <- cur_diff_citestResults
    cur_ord_pag_list[[i]]$curord <- ord
    cur_ord_pag_list[[i]]$scores <- list(ord_symm_diff_score=c(0,0,0),
                                         ord_union_mec_score=c(0,0,0),
                                         mec_score=c(0,0,0))#,
                                         #ord_symm_diff_mse=1,
                                         #mec_mse=1)
    if (i %in% valid_pags_ids) {
      if (!is.null(cur_diff_citestResults) &&
          nrow(cur_diff_citestResults) > 0) {
        indep_probs <- subset(cur_diff_citestResults, type == "indep")$pH0

        dep_probs <- subset(cur_diff_citestResults, type == "dep")$pH1

        unk_probs <- c()

        # symm diff scores:
        ord_symm_diff_score <-
          getProbConjunction(c(dep_probs, indep_probs, unk_probs))
        ord_symm_diff_mse <- getNormalizedSquaredL2DistanceFromCertainty(
          c(dep_probs, indep_probs, unk_probs))
        ord_symm_diff_tuple <- c(ord_symm_diff_score[1],  # upper bound
                                 1 - ord_symm_diff_mse,   # so that the higher, the better
                                 ord_symm_diff_score[2])  # lower bound

        # union scores -- actually uses only CI tests required while
        # computing any of the MECs in the list
        ord_union_mec_score_out <- getStraightforwardPAGScore(
          cur_ord_pag_list[[i]]$amat.pag, unique_ci_relations[, 1:7], m.max = ord)
        ord_union_mec_tuple <- c(ord_union_mec_score_out$score[1], # lowwer bound
                                 1 - ord_union_mec_score_out$mse,  # so that the higher, the better
                                 ord_union_mec_score_out$score[2]) # upper bound

        # mec scores:
        mec_score_tuple <- c(cur_ord_pag_list[[i]]$mec$mec_score[1], # lowwer bound
                             1 - cur_ord_pag_list[[i]]$mec$mec_mse,  # so that the higher, the better
                             cur_ord_pag_list[[i]]$mec$mec_score[2]) # upper bound



        cur_ord_pag_list[[i]]$scores <- list(ord_symm_diff_score = ord_symm_diff_tuple, #ord_symm_diff_score,
                                             ord_union_mec_score = ord_union_mec_tuple,
                                             #ord_symm_diff_mse = ord_symm_diff_mse,
                                             mec_score = mec_score_tuple) #cur_ord_pag_list[[i]]$mec$mec_score,
                                             #mec_mse = cur_ord_pag_list[[i]]$mec$mec_mse)
      } else {
        # mec scores:
        mec_score_tuple <- c(cur_ord_pag_list[[i]]$mec$mec_score[1],
                             1 - cur_ord_pag_list[[i]]$mec$mec_mse, # so that the higher, the better
                             cur_ord_pag_list[[i]]$mec$mec_score[2])

        cur_ord_pag_list[[i]]$scores <- list(ord_symm_diff_score=c(1,1,1),
                                             ord_union_mec_score = c(1,1,1),
                                             #ord_symm_diff_mse=0,
                                             mec_score=mec_score_tuple) #cur_ord_pag_list[[i]]$mec$mec_score,
                                             #mec_mse = cur_ord_pag_list[[i]]$mec$mec_mse)
      }
    }
    cur_ord_pag_list[[i]]$ordPAGs[[as.character(ord)]] <- cur_ord_pag_list[[i]]
  }

  #lapply(cur_ord_pag_list, function(x) { x$scores })

  #lapply(cur_ord_pag_list, function(x) { x$sepsetResults })

    return(list(cur_ord_pag_list=cur_ord_pag_list,
              symmDiffCITests=symmDiffCITests))
}


getTopPagIds <- function(ord_pag_list_score_df,
                         mec_score_df = NULL,
                         ord, sel_top, prob_sel_top,
                         combine_mse = TRUE) {

  top_dc_pag_ids <- NULL

  if (!prob_sel_top) {
    # getting the top "sel_top" pags according to the ord_symm_diff PAG score.
    top_ord_score_pag_ids <- ord_pag_list_score_df[
      which(ord_pag_list_score_df$violation == FALSE &
              ord_pag_list_score_df$index <= sel_top), "pag_list_id"]
    if (!is.null(mec_score_df)) {
      # getting the top "sel_top" pags according to the mec PAG score.
      top_mec_score_pag_ids <- mec_score_df[
        which(mec_score_df$violation == FALSE &
              mec_score_df$index <= sel_top), "pag_list_id"]
    }
  } else if (prob_sel_top >= 1)  {
    top_ord_score_pag_ids <- ord_pag_list_score_df[
      which(ord_pag_list_score_df$violation == FALSE &
              ord_pag_list_score_df$prob_index <= as.numeric(prob_sel_top)), "pag_list_id"]
    if (!is.null(mec_score_df)) {
      top_mec_pag_ids <- mec_score_df[
        which(mec_score_df$violation == FALSE &
              mec_score_df$prob_index <= as.numeric(prob_sel_top)), "pag_list_id"]
    }
  } else {
    stop("Use a numeric value for prob_sel_top to select PAGs based on their score interval.")
  }

  if (combine_mse) {
    agg_scores_up <- c("agg_symmdiff_up_1mse")
    if (!is.null(mec_score_df)) {
      agg_scores_up <- c(agg_scores_up, "agg_mec_up_1mse")
    }
  } else {
    agg_scores_up <- c("agg_symmdiff_up", "agg_symmdiff_1mse")
    if (!is.null(mec_score_df)) {
      agg_scores_up <- c(agg_scores_up, "agg_mec_up", "agg_mec_1mse")
    }
  }

  min_union_score_up <- min(subset(ord_pag_list_score_df, pag_list_id %in% top_ord_score_pag_ids)[, paste0("ord", ord, "_union_score_up")])
  min_mec_score_up <-  min(subset(mec_score_df, pag_list_id %in% top_mec_score_pag_ids)[, paste0("ord", ord, "_mec_score_up")])

  top_union_score_up_ids <- ord_pag_list_score_df[which(ord_pag_list_score_df[, paste0("ord", ord, "_union_score_up")]
                          >= min_union_score_up), "pag_list_id"]
  top_mec_score_up_ids <- mec_score_df[which(mec_score_df[, paste0("ord", ord, "_mec_score_up")]
                                                        >= min_mec_score_up), "pag_list_id"]

  top_ord_pag_ids <- intersect(top_union_score_up_ids, top_mec_score_up_ids)
  top_dc_pag_ids <- unique(c(top_ord_pag_ids, top_mec_score_pag_ids, top_ord_score_pag_ids))

  min_score_considered <- min(subset(ord_pag_list_score_df, pag_list_id %in% top_dc_pag_ids)[, "agg_symmdiff_up_1mse"])
  top_dc_pag_ids <- unique(c(top_dc_pag_ids,
                           ord_pag_list_score_df[which(ord_pag_list_score_df$agg_symmdiff_up_1mse >= min_score_considered), "pag_list_id"]))

  top_dc_pag_ids <- subset(ord_pag_list_score_df,
                           pag_list_id %in% top_dc_pag_ids &
                           violations == FALSE &
                           duplicated == FALSE)$pag_list_id
  return(top_dc_pag_ids)
}




#' @export scorePAGFromCITests
scorePAGFromCITests <- function(amat.pag, ord=Inf, citestResults, verbose=FALSE) {
  a_citestResults <- getAnnotatedCITestResults(amat.pag, ord, citestResults,
                                               verbose=verbose)

  score <- getProbConjunction(c(subset(a_citestResults, type == "indep")$pH0,
                                subset(a_citestResults, type == "dep")$pH1))
  return(list(score=score, a_citestResults=a_citestResults))
}


#' @export getSxyIndepScore
getSxyIndepScore <- function(x, y, Sxy_list, labels,
                             citestResults, indepTest, suffStat,
                             test_only_max=FALSE, # similar to OR / disjunction
                             verbose=FALSE, allowNewTests=TRUE) {
  # getting separation scores
  indep_pH0s <- c()
  for (Sxy in Sxy_list) {
    citest_out <- getCITestResults(x, y, Sxy, citestResults, indepTest,
                                   suffStat, verbose, allowNewTests)
    citestResults <- citest_out$citestResults
    indep_pH0s <- c(indep_pH0s, citest_out$pH0)
  }
  if (verbose && length(indep_pH0s) > 1) {
    cat("There are more than one minimal separator.\n")
    if (test_only_max) {
      cat("Using only the one with the highest score. \n")
    }
  }

  if (test_only_max) {
    indep_pH0s <- max(indep_pH0s)
  }
  indep_score_int <- getIndepScore(h0_probs = indep_pH0s, ret_int=TRUE)

  return(list(indep_pH0=indep_pH0s,
              indep_score_int=indep_score_int,
              citestResults=citestResults))
}


#' @export getSxyMinimalityScore
getSxyMinimalityScore <- function(x, y, Sxy_list, labels, citestResults,
                                  indepTest, suffStat, verbose=FALSE) {
  min_score <- c(1, 1)

  resultsxy <- subset(citestResults, X == x & Y == y)
  resultsyx <- subset(citestResults, X == y & Y == x)
  cur_citests <- rbind(resultsxy, resultsyx)

  cur_subsets_citests <- c()
  for (Sxy in Sxy_list) {
    #cur_citest <- subset(citestResults, X == x & Y == y & S==getSepString(Sxy))
    cur_subsets_citests <- rbind(cur_subsets_citests,
                                 cur_citests[ sapply(cur_citests$S, isProperSubset,
                                               setStr=Sxy, all_subsets=TRUE), ])
  }
  cur_subsets_citests <- cur_subsets_citests[which(!duplicated(cur_subsets_citests)), ]

  min_score_list <- list()
  h1_probs <- cur_subsets_citests$pH1
  if (length(h1_probs) > 0) {
    min_score_list[[length(min_score_list) + 1]] <- getDepScore(h1_probs)
  }

  if (length(min_score_list) > 0) {
    min_score <- getProbConjunctionFromIntervals(min_score_list)
  }
  return(list(min_score=min_score, min_score_list=min_score_list, dep_citests=cur_subsets_citests))
}

#' @importFrom jsonlite toJSON
addTripletListScores <- function(tripletList, sepset, labels, citestResults,
                                 indepTest, suffStat, digits=10, verbose=FALSE,
                                 allowNewTests=TRUE) {
  dep_citests <- data.frame()
  if (!is.null(tripletList)) {
    for (i in 1:nrow(tripletList)) {
      x0 <- tripletList[i, "X0"]
      y0 <- tripletList[i, "Y0"]

      #######################################
      # Getting separation stats and scores #
      #######################################

      # Note: these should be already accounted for in the skeleton score

      SepX0Y0 <- sepset[[x0]][[y0]] # can be a list of minimal sepsets
      if (!is.null(SepX0Y0) && !is.list(SepX0Y0)) {
        SepX0Y0 <- list(SepX0Y0)
      }
      tripletList[i, "SepX0Y0"] <- toJSON(lapply(SepX0Y0, getSepString))
      indep_stats_out <- getSxyIndepScore(x0, y0, SepX0Y0, labels, citestResults,
                                          indepTest, suffStat, verbose)
      citestResults <- indep_stats_out$citestResults

      # can be a list of the pH0 associated to each minimal sepset Sxy
      tripletList[i, "indep_pH0"] <- toJSON(indep_stats_out$indep_pH0,
                                            digits = digits)
      # the interval of the conjunction of all pH0s
      tripletList[i, "indep_score"] <- toJSON(indep_stats_out$indep_score_int,
                                              digits = digits)

      ##############################################
      # Getting each connection's stats and scores #
      ##############################################

      # These should be already accounted for in the case of non-colliders triplets,
      # but may include new ci info in the case of collider triplets!

      # sepset that leaves the path with the triplet open
      Sdep <- getSepVector(tripletList[i, "dep_sepset"])
      npaths <- tripletList[i, "nconnpaths"]

      citest_out <- getCITestResults(x0, y0, Sdep, citestResults, indepTest,
                                     suffStat, verbose, allowNewTests)
      citestResults <- citest_out$citestResults
      tripletList[i, "dep_pH0"] <- citest_out$pH0 # only one numeric value
      dep_score_int <- getDepScore(h1_probs = citest_out$pH1, ret_int=TRUE)
      tripletList[i, "dep_score"] <- toJSON(dep_score_int, digits = digits)

      # These are surely in citestResults
      cur_ciinfo <- subset(citestResults,  X == x0 & Y == y0 & S == getSepString(Sdep))
      cur_ciinfo <- rbind(cur_ciinfo,
                          subset(citestResults,  X == y0 & Y == x0 & S == getSepString(Sdep)))
      cur_ciinfo <- cur_ciinfo[1, ]

      # adding current test to dep_citests
      dep_citests <- rbind(dep_citests,
                         cbind(cur_ciinfo, type = "dep", structure="triplet"))

    }
  }
  return(list(tripletList=tripletList, citestResults=citestResults,
              dep_citests=dep_citests))
}

# This corresponds to the score for all adjacency tests
scoreSkel <- function(sepset, ord, citestResults, indepTest,
                      suffStat, verbose, allowNewTests=TRUE) {
  nvars <- length(sepset)
  tested_ci <- data.frame()
  for (i in 1:(nvars-1)) {
    for (j in (i+1):nvars) {
      Sxy_list <- sepset[[i]][[j]]
      if (is.null(Sxy_list)) {
        # there is an edge, so all tests up to order size.max should indicate dep

        Sxy_list <- getSubsets(aset=1:nvars,
                               only_proper = FALSE, size.max = ord)
        Sxy_list <- Sxy_list[sapply(Sxy_list, function(x) { !(any(c(i,j) %in% x)) })]

        for (Sxy in Sxy_list) {
          cur_citest_out <- getCITestResults(i, j, Sxy, citestResults, indepTest,
                                             suffStat, verbose)
          # adding current test to citestResults
          citestResults <- cur_citest_out$citestResults
          cur_ciinfo <- subset(citestResults,  X == i & Y == j & S == getSepString(Sxy))

          # adding current test to tested_ci
          tested_ci <- rbind(tested_ci,
                             cbind(cur_ciinfo, type = "dep", structure="edge"))
        }
      } else {
        if (!is.list(Sxy_list)) {
          Sxy_list <- list(Sxy_list)
        }
        # there is one or more minimal separating sets
        for (Sxy in Sxy_list) {
          # adding ci info associated to the minimal independence
          cur_citest_out <- getCITestResults(i, j, Sxy, citestResults, indepTest,
                                             suffStat, verbose)
          # adding current test to citestResults
          citestResults <- cur_citest_out$citestResults
          cur_ciinfo <- subset(citestResults,  X == i & Y == j & S == getSepString(Sxy))

          # adding current test to tested_ci
          tested_ci <- rbind(tested_ci,
                             cbind(cur_ciinfo, type = "indep", structure="separ"))

          # adding ci info associated to dep given proper subsets of the minimal separator
          pSxy_list <- getSubsets(aset=Sxy, only_proper = TRUE)
          for (pSxy in pSxy_list) {
            cur_citest_out <- getCITestResults(i, j, pSxy, citestResults, indepTest,
                                               suffStat, verbose)
            # adding current test to citestResults
            citestResults <- cur_citest_out$citestResults
            cur_ciinfo <- subset(citestResults,  X == i & Y == j & S == getSepString(pSxy))

            # adding current test to tested_ci
            tested_ci <- rbind(tested_ci,
                               cbind(cur_ciinfo, type = "dep", structure="min"))
          }
        }
      }
    }
  }
  tested_ci <- tested_ci[order(tested_ci$ord), ]
  citestResults <- citestResults[order(citestResults$ord), ]

  skel_score <- getProbConjunction(
    c(subset(tested_ci, type == "indep")$pH0,
      subset(tested_ci, type == "dep")$pH1))

  return(list(skel_citests=tested_ci, skel_score=skel_score,
              citestResults=citestResults))
}

#' @export scoreMEC
scoreMEC <- function(mec, sepset, max.ord, citestResults, indepTest,
                     suffStat, verbose, allowNewTests=TRUE) {
  scored_mec <- mec
  labels <- colnames(mec$skel)

  out_skel <- scoreSkel(sepset, max.ord, citestResults, indepTest,
                        suffStat, verbose, allowNewTests=allowNewTests)
  scored_mec$skel_citests <- out_skel$skel_citests
  scored_mec$skel_score <- out_skel$skel_score
  citestResults <- out_skel$citestResults

  scored_mec$ck_dep_citests <- c()
  if (!is.null(scored_mec$CK)) {
    out_ck <- addTripletListScores(scored_mec$CK, sepset, labels, citestResults,
                                   indepTest, suffStat,
                                   verbose=verbose, allowNewTests=allowNewTests)
    out_ck$dep_citests[, "structure"] <- "ck"
    citestResults <- out_ck$citestResults
    scored_mec$CK <- out_ck$tripletList
    scored_mec$ck_dep_citests <- out_ck$dep_citests
  }

  scored_mec$nck_dep_citests <- c()
  if (!is.null(scored_mec$NCK)) {
    out_nck <- addTripletListScores(scored_mec$NCK, sepset, labels, citestResults,
                                    indepTest, suffStat,
                                    verbose=verbose, allowNewTests=allowNewTests)
    out_nck$dep_citests[, "structure"] <- "nck"
    scored_mec$NCK <- out_nck$tripletList
    citestResults <- out_nck$citestResults
    scored_mec$nck_dep_citests <- out_nck$dep_citests
  }

  scored_mec$all_citests <- rbind(
    scored_mec$skel_citests,
    scored_mec$ck_dep_citests,
    scored_mec$nck_dep_citests
  )

  scored_mec$all_citests <- scored_mec$all_citests[!duplicated(
    scored_mec$all_citests[, -ncol(scored_mec$all_citests)]) , ]
  scored_mec$all_citests <- scored_mec$all_citests[order(scored_mec$all_citests$ord), ]

  mec_probs <- c(subset(scored_mec$all_citests, type == "indep")$pH0,
                 subset(scored_mec$all_citests, type == "dep")$pH1)

  scored_mec$mec_score <- getProbConjunction(mec_probs)

  scored_mec$mec_mse <- getNormalizedSquaredL2DistanceFromCertainty(mec_probs)

  return(list(scored_mec=scored_mec,
              citestResults=citestResults))
}


# Returns probability interval that all dependencies implied by the separation
# between x0 and y0 are true.
#' @importFrom jsonlite fromJSON
#' @export getTripletDepScores
getTripletDepScores <- function(scored_mec, x0, y0) {
  dep_citests <- data.frame()
  dep_score_list <- list()
  dep_score <- c(NA, NA)
  if (!is.null(scored_mec$CK)) {
    mec_ck <- subset(scored_mec$CK, X0 == x0 & Y0 == y0)
    if (dim(mec_ck)[1] > 0) {
      dep_score_list <- c(dep_score_list,
                          lapply(mec_ck$dep_score, fromJSON))
      dep_citests <- rbind(dep_citests, c(mec_ck[, c("X0", "Y0", "dep_sepset")]))
    }
  }
  if (!is.null(scored_mec$NCK)) {
    mec_nck <- subset(scored_mec$NCK, X0 == x0 & Y0 == y0)
    if (dim(mec_nck)[1] > 0) {
      dep_score_list <- c(dep_score_list,
                          lapply(mec_nck$dep_score, fromJSON))
      dep_citests <- rbind(dep_citests, c(mec_nck[, c("X0", "Y0", "dep_sepset")]))
    }
  }
  if (length(dep_score_list) > 0) {
    dep_score <- getProbConjunctionFromIntervals(dep_score_list)
  }
  return(list(dep_score=dep_score, dep_score_list=dep_score_list, dep_citests=dep_citests))
}


