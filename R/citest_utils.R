getCITestResultsHelper <- function(x, y, Sxy, citestResults) {
  X = Y = S = NULL

  SxyStr <- getSepString(sort(Sxy))
  pvalueXY.S <- pH0XY.S <- pH1XY.S <- NULL
  resultsxys <- subset(citestResults, X == x & Y == y & S == SxyStr)
  resultsyxs <- subset(citestResults, X == y & Y == x & S == SxyStr)
  resultsxys <- rbind(resultsxys, resultsyxs)
  if (!is.null(resultsxys) && nrow(resultsxys) > 0) {
    # the test should be the symmetric for X,Y|S and Y,X|S
    pvalueXY.S <- resultsxys[1, "pvalue"]
    pH0XY.S <- resultsxys[1, "pH0"]
    pH1XY.S <- resultsxys[1, "pH1"]
  }

  if (is.null(pvalueXY.S) || is.null(pH1XY.S) || is.null(pH0XY.S)) {
    return(NULL)
  } else {
    return(list(pvalueXY.S=pvalueXY.S, pH0XY.S=pH0XY.S, pH1XY.S=pH1XY.S))
  }
}

# First, tries to get results from citestResults, which are those
# already computed and saved by the algorithms.
# If those are not available yet, tries to get from suffStat$citestResults
# if this is provided. Otherwise, computes a new test.
#' @importFrom FCI.Utils getSepString
#' @export getCITestResults
getCITestResults <- function(x, y, Sxy, citestResults, indepTest, suffStat,
                             NAdelete, verbose=FALSE, allowNewTests=TRUE) {
  pvalue = pH0 = pH1 = NULL

  sortedxy <- sort(c(x,y))
  x <- sortedxy[1]
  y <- sortedxy[2]
  out <- getCITestResultsHelper(x, y, Sxy, citestResults)
  SxyStr <- getSepString(sort(Sxy))
  curResults <-  data.frame("ord"= length(Sxy),
                            "X"=x, "Y"=y,
                            "S"=SxyStr)

  if (!is.null(out)) {
    pvalue = out$pvalueXY.S
    pH0 = out$pH0XY.S
    pH1 = out$pH1XY.S
  } else if (!is.null(suffStat$citestResults)) {
    # getting results from suffStat$citestResults if it is provided.
    out <- getCITestResultsHelper(x, y, Sxy, suffStat$citestResults)
    if (!is.null(out)) {
      pvalue = out$pvalueXY.S
      pH0 = out$pH0XY.S
      pH1 = out$pH1XY.S
      curResults <- cbind(curResults, "pvalue"=pvalue, "pH0"=pH0, "pH1"=pH1)
      citestResults <- rbind(citestResults, curResults)
    }
  }

  if (is.null(pvalue) || is.null(pH0) || is.null(pH1)) {
    # performing tests only if results have not been provided.
    if (allowNewTests) {
      if (verbose) cat("Performing new test...")

      testXY.S <- indepTest(x, y, Sxy, suffStat)
      pvalue <- testXY.S
      if (is.na(pvalue))
        pvalue <- as.numeric(NAdelete)

      if (is.null(suffStat$n)) {
        stop("Please provide the sample size \'n\' in the suffStat list.")
      }
      probsXY.S <- pvalue2probs(pvalue, n=suffStat$n)
      pH0 <- probsXY.S$pH0
      pH1 <- probsXY.S$pH1

      curResults <- cbind(curResults, "pvalue"=pvalue,
                        "pH0"=pH0, "pH1"=pH1)

      citestResults <- rbind(citestResults, curResults)
    } else {
      errmessage <- paste0("Err: Test for ", x,",", y, "|", SxyStr, " is not available.")
      fileConn <- file(paste0("log.txt"), open = "a")
      writeLines(errmessage, fileConn)
      writeLines(paste0(" Labels=", paste0(colnames(suffStat$dataset), collapse=";")), fileConn)
      close(fileConn)

      stop(errmessage)
    }
  }

  if (verbose) {
    if (!is.null(pvalue) && !is.null(pH0) && !is.null(pH1)) {
      cat("\n stats for ",x,",", y, "|", SxyStr, ":  pvalue =", pvalue,
        "; pH0 =",  pH0, "; pH1 =", pH1, "\n")
    } else {
      cat("\n Cannot compute stats for ",x,",", y, "|", SxyStr, "\n")
    }
  }

  return(list(pvalue=pvalue, pH0=pH0, pH1=pH1, citestResults = citestResults))
}
