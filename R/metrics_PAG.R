#source("mec_utils.R")

#' @export getMECFaithfulnessScores
getMECFaithfulnessScores <- function(apag, suffStat,
                                     alphas = c(0.01, 0.05), ord=NULL) {
  mec_score <- c(NA, NA)
  faithf_mec_dep_counts <- c(NA, NA, NA)
  faithf_mec_indep_counts <- c(NA, NA, NA)

  mec <- getMEC(apag, scored = TRUE, max.ord = ord, citestResults = suffStat$citestResults)
  if (!is.null(mec)) {
    mec_faithfulness <- getFaithfulnessDegree(apag, mec$mec$all_citests)
    mec_ci <- mec_faithfulness$f_citestResults
    mec_dep_pvalues <- subset(mec_ci, type == "dep")$pvalue
    mec_indep_pvalues <- subset(mec_ci, type == "indep")$pvalue

    faithf_mec_dep_counts <- c(colSums( matrix(unlist(
      lapply(alphas, function(x) { mec_dep_pvalues < x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
      length(mec_dep_pvalues))

    faithf_mec_indep_counts <- c(colSums( matrix(unlist(
      lapply(alphas, function(x) { mec_indep_pvalues > x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
      length(mec_indep_pvalues))

    mec_score <- getProbConjunction(mec_faithfulness$probs)
  }
  return(list(mec_score=mec_score, faithf_mec_dep_counts=faithf_mec_dep_counts,
              faithf_mec_indep_counts=faithf_mec_indep_counts))
}

#' @export getStraightforwardPAGScore2
getStraightforwardPAGScore2 <- function(apag, suffStat, ord=NULL) {
  # this checks the agreement with the complete citestResults table
  faithfulnessDegree <- getFaithfulnessDegree(apag, suffStat$citestResults)
  l2d <- getL2DistanceFromCertainty(faithfulnessDegree$probs)
  l2d2 <- getSquaredL2DistanceFromCertainty(faithfulnessDegree$probs)
  mse <- getNormalizedSquaredL2DistanceFromCertainty(faithfulnessDegree$probs)
  frechet_int <- getProbConjunction(faithfulnessDegree$probs)
  return(list(mse=mse, l2d=l2d, l2d2=l2d2, frechetLB=frechet_int[1], frechetUB=frechet_int[2]))
}

#' @export getMECTargetedPAGScore
getMECTargetedPAGScore <- function(apag, suffStat, ord=NULL,
                                   allowNewTests = FALSE,
                                   indepTest = NULL,
                                   verbose = FALSE) {
  mec <- getMEC(apag, scored = TRUE, max.ord = ord,
                citestResults = suffStat$citestResults,
                allowNewTests = allowNewTests, verbose=verbose,
                indepTest = indepTest, suffStat = suffStat)
  citestResults <- mec$citestResults
  faithfulnessDegree <- getFaithfulnessDegree(apag, mec$mec$all_citests)

  l2d <- getL2DistanceFromCertainty(faithfulnessDegree$probs)
  l2d2 <- getSquaredL2DistanceFromCertainty(faithfulnessDegree$probs)
  mse <- getNormalizedSquaredL2DistanceFromCertainty(faithfulnessDegree$probs)
  frechet_int <- getProbConjunction(faithfulnessDegree$probs)
  return(list(mse=mse, l2d=l2d, l2d2=l2d2, frechetLB=frechet_int[1], frechetUB=frechet_int[2],
              citestResults=citestResults))
}


#' @export getFaithfulnessScores
getFaithfulnessScores <- function(apag, suffStat,
                                  alphas = c(0.01, 0.05), ord=NULL) {
  if (is.null(ord)) {
    ord <- ncol(apag) - 2
  }

  faithf_dep_counts <- c(NA, NA, NA)
  faithf_indep_counts <- c(NA, NA, NA)
  faithf_adj_dep_counts <- c(NA, NA, NA)
  faithf_skel_dep_counts <- c(NA, NA, NA)
  faithf_skel_indep_counts <- c(NA, NA, NA)
  faithf_mec_dep_counts <- c(NA, NA, NA)
  faithf_mec_indep_counts <- c(NA, NA, NA)

  # this tests the agreement with the complete citestResults table
  faithfulnessDegree <- getFaithfulnessDegree(apag, suffStat$citestResults)
  trivial_score <- getProbConjunction(faithfulnessDegree$probs)

  dep_pvalues <- subset(faithfulnessDegree$f_citestResults, type == "dep")$pvalue
  indep_pvalues <- subset(faithfulnessDegree$f_citestResults, type == "indep")$pvalue
  faithf_dep_counts <- c(colSums( matrix(unlist(
    lapply(alphas, function(x) { dep_pvalues < x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
    length(dep_pvalues))

  faithf_indep_counts <- c(colSums( matrix(unlist(
    lapply(alphas, function(x) { indep_pvalues > x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
    length(indep_pvalues))

  # this tests the agreement with the adj-required citest results

  adj_score <- c(NA, NA)
  adj_faithfulnessDegree <- getAdjFaithfulnessDegree(
    apag, suffStat$citestResults, test_dep_min = FALSE, ord)
  if (!is.null(adj_faithfulnessDegree)) {
    adj_ci <- adj_faithfulnessDegree$tested_ci
    adj_dep_pvalues <- subset(adj_ci, type == "dep")$pvalue

    faithf_adj_dep_counts <- c(colSums( matrix(unlist(
      lapply(alphas, function(x) { adj_dep_pvalues < x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
      length(adj_dep_pvalues))

    adj_score <- getProbConjunction(c(
      subset(adj_ci, type == "indep", select = pH0, drop=TRUE), # there is never indep
      subset(adj_ci, type == "dep", select = pH1, drop=TRUE)))
  }

  skel_score <- c(NA, NA)
  skel_faithfulnessDegree <- getAdjFaithfulnessDegree(
    apag, suffStat$citestResults, test_dep_min = TRUE, ord)
  if (!is.null(skel_faithfulnessDegree)) {
    skel_ci <- skel_faithfulnessDegree$tested_ci
    skel_dep_pvalues <- subset(skel_ci, type == "dep")$pvalue
    skel_indep_pvalues <- subset(skel_ci, type == "indep")$pvalue

    faithf_skel_dep_counts <- c(colSums( matrix(unlist(
      lapply(alphas, function(x) { skel_dep_pvalues < x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
      length(skel_dep_pvalues))

    faithf_skel_indep_counts <- c(colSums( matrix(unlist(
      lapply(alphas, function(x) { skel_indep_pvalues > x})), ncol = 2, byrow=FALSE), na.rm=TRUE),
      length(skel_indep_pvalues))

    skel_score <- getProbConjunction(c(
      subset(skel_ci, type == "indep", select = pH0, drop=TRUE), # there is never indep
      subset(skel_ci, type == "dep", select = pH1, drop=TRUE)))
  }

  mec_score_out <- getMECFaithfulnessScores(apag, suffStat, alphas = alphas, ord=ord)
  mec_score <- mec_score_out$mec_score
  faithf_mec_dep_counts <- mec_score_out$faithf_mec_dep_counts
  faithf_mec_indep_counts <- mec_score_out$faithf_mec_indep_counts

  return(list(trivial_score=trivial_score,
              adj_score=adj_score,
              mec_score=mec_score,
              faithf_dep_counts=faithf_dep_counts,
              faithf_indep_counts=faithf_indep_counts,
              faithf_adj_dep_counts=faithf_adj_dep_counts,
              faithf_skel_dep_counts=faithf_skel_dep_counts,
              faithf_skel_indep_counts=faithf_skel_indep_counts,
              faithf_mec_dep_counts=faithf_mec_dep_counts,
              faithf_mec_indep_counts=faithf_mec_indep_counts))
}

# Only set checkViolations = FALSE if est.amat.pag is valid and has no violations
#' @export getMetrics
getMetrics <- function(true.amat.pag, est.amat.pag, est.sepset=NULL, dat=NULL,
                       conservative=FALSE, checkViolations = TRUE) {
  fdr=NA; fomr=NA; shd=NA
  if (!is.null(est.amat.pag)) {
    posneg <- getPAGPosNegMetrics(true.amat.pag, est.amat.pag)
    fdr <- posneg$false_discovery_rate
    fomr <- posneg$false_omission_rate
    shd = shd_PAG(true.amat.pag, est.amat.pag)
  }

  checkPAGValidity <- TRUE
  if (!checkViolations) {
    checkPAGValidity <- FALSE
  }

  bic <- NA
  hasViol <- msep_violations <- valid_PAG <- NA
  if (checkViolations && !is.null(est.sepset)) {
    violations <- hasViolation(est.amat.pag, sepset=est.sepset,
                               conservative=conservative,
                               log=TRUE, verbose=FALSE)
    valid_PAG <- violations$log$validPAG
    hasViol <- violations$out
    if (!is.null(violations$log$`m-separations`$msep) &&
        !is.null(violations$log$`def_m-connections_min`$mconn == FALSE)) {
      msep_violations <- length(which(violations$log$`m-separations`$msep == FALSE)) +
        length(which(violations$log$`def_m-connections_min`$mconn == FALSE))
    } else {
      msep_violations <- NA
    }
  } else if (checkPAGValidity) {
    valid_PAG <- isValidPAG(est.amat.pag, conservative = conservative,
                            knowledge = FALSE, verbose=FALSE)
  }

  if ((!checkPAGValidity || valid_PAG) && !is.null(dat)) {
    bic <- getBIC(est.amat.pag, dat, type="BIC")$bic
  }

  return(data.frame(fdr=fdr, fomr=fomr, shd=shd, bic=bic,
                    valid=valid_PAG, viol=hasViol,
                    msepViol=msep_violations))
}

