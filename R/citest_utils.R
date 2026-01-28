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

    if (is.na(pvalueXY.S) || is.na(pH1XY.S) || is.na(pH0XY.S)) {
        pH0XY.S <- 0.5
        pH1XY.S <- 0.5
    }
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
                             verbose=FALSE, allowNewTests=TRUE) {
  pvalue = pH0 = pH1 = NULL
  chi2stat = df = NULL # default will be used for LR GLM-based CI Test

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

      if (suffStat$retall) {
        pvalue <- testXY.S$p
        if (suffStat$method == "nnGCM") {
          pvalue <- testXY.S$ret$pvalue
          chi2stat <- testXY.S$ret$chi2stat
          df <- testXY.S$ret$df
        }
      } else {
        pvalue <- testXY.S
      }

      if (is.null(suffStat$n)) {
        if (!is.null(suffStat$dataset)) {
          suffStat$n <- nrow(suffStat$dataset)
        } else {
          stop("Please provide the sample size \'n\' in the suffStat list.")
        }
      }

      if (is.na(pvalue)) {
        pH0 <- 0.5
        pH1 <- 0.5
      } else {
        probsXY.S <- pvalue2probs(pvalue, n=suffStat$n, eff_size=suffStat$eff_size,
                              chiSqStat = chi2stat, df=df)
        pH0 <- probsXY.S$pH0
        pH1 <- probsXY.S$pH1
      }

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
