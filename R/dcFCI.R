getOrdPotSepsets <- function(pagAdjM, ord, alpha, amec,
                             fixedEdges, doneEdges, citestResults,
                             indepTest, suffStat, verbose,
                             pH0ThreshMin = 0.5,
                             allowNewTests=TRUE) {
  bestSepResults <- c()

  adjM <- pagAdjM > 0
  # indices of the pairs of connected variables
  ind <- as.data.frame(which(lower.tri(adjM) & adjM, arr.ind = TRUE))
  rownames(ind) <- NULL

  # iterates over all edges
  remEdge_i <- 1
  #TODO: make this more efficient / parallelize
  while (remEdge_i <= nrow(ind)) {
    curEdgeResults <- c()
    x <- as.numeric(ind[remEdge_i, 1])
    y <- as.numeric(ind[remEdge_i, 2])
    # continues only if there is an edge between the pair of variables
    # and such an edge is not in the set of edges that cannot be removed
    if (pagAdjM[y, x] && !fixedEdges[y, x] && !doneEdges[y,x]) {
      potSepsOut <- getOrdPotSepsetsXY(x, y, ord, pagAdjM, verbose)
      if (potSepsOut$done) {
        doneEdges[x,y] <- doneEdges[y,x] <- TRUE
        if (verbose) {
          cat(paste0("Edge ", x, "-", y, " is done.\n"))
        }
      }
      potSeps <- potSepsOut$potsepsets
      n_potSeps <- length(potSeps)

      if (!is.null(potSeps)) {
        if (verbose)
          print(paste0("potSeps for x=", x, " y=", y, " ord=", ord, " is {",
                       paste0(potSeps, collapse="; "), "}"))

        for (Si in potSeps) {
          citest_out <- getCITestResults(x, y, Si, citestResults, indepTest,
                                         suffStat, verbose, allowNewTests)
          pvalue <- citest_out$pvalue
          pH0 <- citest_out$pH0
          pH1 <- citest_out$pH1
          citestResults <- citest_out$citestResults

          if (is.na(pvalue) || pvalue > alpha || pH0 >= pH0ThreshMin) { # PAGs with and without the edge will be considered
            if (verbose) {
              cat(paste0("Adding ", getSepString(Si), " to potsepsets of ", x, " and ", y, "\n"))
            }

            curEdgeResults <- rbind(curEdgeResults,
                                    data.frame("ord"= length(Si),
                                               "X"=x, "Y"=y,
                                               "S"=getSepString(Si),
                                               "pvalue"=pvalue,
                                               "pH0"=pH0,
                                               "pH1"=pH1))

          }
        }
        bestSepResults <- rbind(bestSepResults,
                                curEdgeResults)
      } else {
        if (verbose)
          print(paste0("potSeps for x=", x, " y=", y, " ord=", ord, " is NULL"))
      }
    }
    remEdge_i <- remEdge_i + 1
  }

  # We now sort using pH0
  if (!is.null(bestSepResults) && nrow(bestSepResults) > 0) {
    bestSepResults <- bestSepResults[order(bestSepResults$pH0, decreasing = T),]
  }


  return(list(sepsetResults=bestSepResults, citestResults=citestResults,
              doneEdges=doneEdges))
}


# This gets the potential sepsets of a pair of variables (X,Y) by looking
# at all shortest possibly m-connecting paths in a PAG constructed
# without the edge X-Y and only with the faithful colliders
getOrdPotSepsetsXY <- function(x, y, ord, pagAdjM, verbose=FALSE) {
  labels <- colnames(pagAdjM)
  p <- dim(pagAdjM)[1]

  mutPagAdjM <- pagAdjM
  mutPagAdjM[x,y] <- mutPagAdjM[y,x] <- 0
  potSeps <- list()

  done <- ord > (p-2)
  if (!done) {
    potSeps[[1]] <- vector("numeric", 0)
    curord <- 1
    while (!done && curord <= ord) {
      potssid = 1
      while (potssid <= length(potSeps) && length(potSeps[[potssid]]) < curord) {
        curpotsepids <- potSeps[[potssid]]
        if (length(curpotsepids) == 0) {
          potSepNames = c()
        } else {
          potSepNames <- labels[curpotsepids]
        }
        mconnpaths <- getMConnPaths(mutPagAdjM, labels[x], labels[y], potSepNames,
                                    definite = FALSE,
                                    type = "all", maxLength = curord+2,
                                    verbose=verbose)
        potSeps <- potSeps[-potssid]

        if (length(mconnpaths) == 0 && length(potSeps) == 0) {
          if (verbose) {
            cat(paste0("It is not possible to find a sepset for {", labels[x], ",", labels[y],
                       "} of order ", curord, " (<=", ord, ") or higher\n"))
          }
          done <- TRUE
          break
        }
        updated <- FALSE
        for (mconnpath in mconnpaths) {
          varid <- 2
          while (varid < min(length(mconnpath), ord+2)) {
            var <- which(labels == mconnpath[varid])
            trplt <- mconnpath[(varid-1):(varid+1)]

            if (verbose) {
              cat("Checking triplet", paste0(trplt, collapse="-"), "\n")
            }

            if (!(mconnpath[varid] %in% potSepNames)) {
              if (!isCollider(mutPagAdjM, trplt))  {
                potSeps[[length(potSeps)+1]] <- sort(c(curpotsepids, var))
                updated <- TRUE
                break
              }
            }
            varid = varid + 1
          }
        }
        if (!updated) {
          potssid = potssid + 1
        }
      }
      curord <- curord + 1
    }
  }

  finalPotSeps <- unique(potSeps)

  if (ord == 0 || (length(finalPotSeps) > 0 && length(finalPotSeps[[1]]) == ord)) {
    potsepsets = finalPotSeps
  } else {
    potsepsets = NULL
  }
  return(list(potsepsets = potsepsets, done = done))
}

#' @importFrom jsonlite toJSON fromJSON
getPAGListScores <- function(pag_List, ord_pag_score="ord_symm_diff_score",
                             score_len = 3) {
  violations <- sapply(pag_List, function (x) {x$violations})

  ord_pag_scores <- data.frame()
  for (pag in pag_List) {
    cur_ord_pag_scores <- as.character(toJSON(
      lapply(pag$ordPAGs, function(x) {
        if (!is.null(x$scores)) {
          x$scores[[ord_pag_score]]
        } else {
          rep(0, score_len)
        }
      }), digits=10))
    ord_pag_scores <- rbind.data.frame(ord_pag_scores, data.frame(cur_ord_pag_scores))
  }

  scores_df = cbind.data.frame(violation=violations, ord_pag_scores)
  colnames(scores_df) <- c("violation", "cur_ord_pag_scores")

  return(scores_df)
}


updateSkelScore <- function(apag, ord, citestResults, indepTest,
                            suffStat, verbose,
                            allowNewTests) {
  out_skel <- scoreSkel(apag$sepset, ord,
                        citestResults, indepTest,
                        suffStat, verbose > 1,
                        allowNewTests=allowNewTests)
  apag$mec$skel_citests <- out_skel$skel_citests
  apag$mec$skel_score <- out_skel$skel_score
  apag$mec$all_citests <- unique(rbind(
    apag$mec$all_citests, apag$mec$skel_citests))

  mec_probs <- c(subset(apag$mec$all_citests, type == "indep")$pH0,
                 subset(apag$mec$all_citests, type == "dep")$pH1)

  apag$mec$mec_score <- getProbConjunction(mec_probs)
  apag$mec$mec_mse <- getNormalizedSquaredL2DistanceFromCertainty(mec_probs)

  citestResults <- out_skel$citestResults

  return(list(apag=apag, citestResults=citestResults))
}


processPotPAG <- function(id, toProcessPAGList, ord, citestResults, indepTest,
                          suffStat, verbose=FALSE, allowNewTests=TRUE,
                          appender_file = NULL) {

  futile.logger::flog.appender(futile.logger::appender.file(appender_file),
                               name = "dcfci_log")
  futile.logger::flog.info(
    paste("Started processing PAG", id, " / ", length(toProcessPAGList)), name = "dcfci_log")

  potPAG <- toProcessPAGList[[id]]
  ######################
  # Applying FCI rules #
  ######################
  rules <- rep(TRUE, 10)

  pot_amat.pag <- pcalg::udag2pag(pag = potPAG$amat.pag, potPAG$sepset, rules = rules,
                                  unfVect = c(), verbose = verbose > 1, orientCollider = TRUE)
  potPAG$amat.pag <- pot_amat.pag


  ###########################
  # Checking for violations #
  ###########################
  potPAG$violations <- hasViolation(potPAG$amat.pag, potPAG$sepset,
                                    listall = FALSE, verbose=verbose > 1)

  ########################
  # Computing Scored MEC #
  ########################
  if (!potPAG$violations) {
    if (verbose)
      cat("Updating MEC of new candidate PAG...\n")
    out_mec <- getMEC(potPAG$amat.pag, ag.type="pag", scored=TRUE,
                      max.ord = ord, citestResults,
                      indepTest, suffStat,
                      verbose=verbose > 1, allowNewTests=allowNewTests)
    potPAG$mec <- out_mec$mec
    citestResults <- out_mec$citestResults
  } else {
    potPAG$mec <- NULL
  }

  potPAG$done <- TRUE
  futile.logger::flog.info(
    paste("Finished processing PAG", id, " / ", length(toProcessPAGList)), name = "dcfci_log")


  return(list(potPAG=potPAG, citestResults=citestResults))
}


#' @importFrom futile.logger flog.appender appender.file flog.info
#' @importFrom rje powerSet
#' @importFrom jsonlite fromJSON
#' @importFrom future.apply future_lapply
#' @importFrom pcalg udag2pag
#' @export dcFCI
dcFCI <- function(suffStat, indepTest, labels, alpha=0.05,
                  m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                  verbose = TRUE, sel_top = 1, prob_sel_top = FALSE,
                  run_parallel = TRUE, allowNewTests=TRUE,
                  list.max = 500,
                  pH0ThreshMin=0.5, pH0ThreshMax=1,
                  order_by_mse = FALSE,
                  log_folder = file.path(getwd(), "tmp", "logs")) {

  ord_pag_score = "ord_symm_diff_score"

  if (is.null(prob_sel_top)) {
    prob_sel_top = 1
  }

  if (is.null(sel_top)) {
    sel_top = 1
  }

  getToProcessPAGList <- function(curpag, ord, fixedEdges, alpha, citestResults,
                                  indepTest, suffStat, verbose=TRUE,
                                  pH0ThreshMin=0.5, pH0ThreshMax=1,
                                  allowNewTests=TRUE) {

    # inner function to generate the sepsetResultsList in parallel
    checkSepsetResults <- function(ids, sepsetResults, certainSepsetResults) {
      #cat("ids:", ids, "\n")
      curSepsetResults <- sepsetResults[ids, ]
      if (nrow(unique(curSepsetResults[,c(2:3)])) < nrow(curSepsetResults)) {
        # there are "duplicated" pairs, so we skip it.
        return(NULL)
      }
      return(rbind(certainSepsetResults, curSepsetResults))
    }

    # inner function to generate all potPAGs in parallel
    createToProcessPAGList <- function(curSepsetResults, curpag, sepset, ord) {
      # generate a PAG for a given curSepsetResults
      # (an element of the sepsetResults of the list

      potAdjM <- curpag$amat.pag > 0
      potSepset <- sepset
      doneEdges <- curpag$doneEdges
      if (!is.null(curSepsetResults) && nrow(curSepsetResults) > 0) {
        for (i in 1:nrow(curSepsetResults)) {
          x = curSepsetResults[i,"X"]
          y = curSepsetResults[i,"Y"]
          S = getSepVector(curSepsetResults[i, "S"])

          potAdjM[x, y] <- potAdjM[y, x] <- 0
          potSepset[[x]][[y]] <- potSepset[[y]][[x]] <- S
          doneEdges[x, y] <- doneEdges[y, x] <- TRUE
        }
      }

      potPAG <-  list(amat.pag = potAdjM, # to be processed
                      sepset = potSepset, # the used sepset
                      mec = NULL, # to be computed
                      scores = curpag$scores,
                      sepsetResults = rbind(curpag$sepsetResults, curSepsetResults), # used sepset with pvalues
                      doneEdges = doneEdges,
                      curord = ord,
                      violations = NA,
                      done = FALSE,
                      ordPAGs=curpag$ordPAGs)
      return(potPAG)
    }


    pagAdjM <-  curpag$amat.pag
    p <- ncol(pagAdjM)
    doneEdges <- curpag$doneEdges
    sepset <- curpag$sepset
    amec <- curpag$amec

    if (verbose) {
      cat("Getting potsepsets...\n")
    }
    sepsetResultsOut <- getOrdPotSepsets(pagAdjM, ord, alpha, amec,
                                         fixedEdges, doneEdges, citestResults,
                                         indepTest, suffStat, verbose > 1,
                                         pH0ThreshMin=pH0ThreshMin,
                                         allowNewTests=allowNewTests)

    citestResults <- sepsetResultsOut$citestResults
    sepsetResults <- sepsetResultsOut$sepsetResults

    if (verbose) {
      if (!is.null(sepsetResults)) {
        cat("Number of sepsets: ", nrow(sepsetResults), "...\n")
      } else {
        cat("Number of sepsets: ", 0, "...\n")
      }
    }

    # updating the input pag (i.e, curpag) in the output list.
    curpag$doneEdges <- sepsetResultsOut$doneEdges
    curpag$done <- TRUE # its sepset won't change

    # # this creates a pag corresponding to each possible combination of the sepsets
    toProcessPAGList <- list()

    # sepsetResults can still have pairs with more than one minimal sepset
    # we have to make the list of all possible sepsetResults, but each containing at most
    # one minimal sepset by pair. Those are already part of the powerSet. So we simply
    # discard sepsetResults with "duplicated" pairs.
    sepsetResultsList <- list()
    if (!is.null(sepsetResults) && nrow(sepsetResults) > 0) {

      # Note: certainsepsets are not only those with pH0 > pH0ThreshMax, but
      # also a potential sepset. Those that are not potsepset are simply ignored.
      certainSepsetResults <- subset(sepsetResults, pH0 >= pH0ThreshMax)

      # if more than one minimal separator is certain, we just get one of them
      # and ignore the others (which will be evaluated in the score as  correct)
      # TODO: this may lead to violations -- should we not consider pairs
      # with more than one "certain" minimal separator as certain?

      certainPairs <- unique(certainSepsetResults[,c(2:3)])
      if (verbose) {
        cat("certainSepsetResults: \n\n")
        print(certainSepsetResults)
        #cat("Certain pairs: \n\n")
        #print(certainPairs)
      }

      # separated pairs from 'certainPairs' should be removed, regardless of the minimal sep
      for (pair_row in 1:nrow(certainPairs)) {
        rows_to_rm <- which(sepsetResults$X==certainPairs[pair_row,1] &
                              sepsetResults$Y == certainPairs[pair_row,2])
        if (length(rows_to_rm) > 0) {
          sepsetResults <- sepsetResults[-rows_to_rm, ]
        }
      }
      if (verbose) {
        cat("sepsetResults after: \n\n")
        print(sepsetResults)
      }


      if (nrow(sepsetResults) > 0) {
        if (verbose) {
          cat(paste0("Number of uncertain sepsets: ", nrow(sepsetResults), "\n"))
          if (nrow(sepsetResults) > 20) {
            cat("ERROR: Unable to proceed due to excessive uncertainty in independence test results.\n\n")
            print(sepsetResults)
            return(NULL)
          }
        }

        #if (nrow(certainSepsetResults) > 0) {
          psSepsetResults <- powerSet(1:nrow(sepsetResults))  # keeps the empty set
        #} else {
        #  psSepsetResults <- powerSet(1:nrow(sepsetResults))[-1]  # removes the empty set
        #}

        if (verbose) {
          cat(paste0("Length of the powerSet of sepsetResults: ", length(psSepsetResults), "\n"))
        }

        sepsetResultsList <- future_lapply(psSepsetResults, checkSepsetResults,
                                           sepsetResults, certainSepsetResults,
                                           future.seed=TRUE)

        sepsetResultsList <- sepsetResultsList[!sapply(sepsetResultsList, is.null)]
      } else {
        sepsetResultsList <- list(certainSepsetResults)
      }

      if (verbose) {
        cat(paste0("Length of sepsetResultsList: ", length(sepsetResultsList), "\n\n"))
      }

      toProcessPAGList <- future_lapply(sepsetResultsList, createToProcessPAGList,
                                        curpag, sepset, ord, future.seed=TRUE)
    }

    return(list(curpag=curpag,
                toProcessPAGList=toProcessPAGList,
                citestResults=citestResults))
  }

  elapsed_time <- system.time({
    rules = rep(TRUE, 10)
    p <- length(labels)

    ## creates the initial skeleton based on fixedGaps and fixedEdges #
    ###################################################################

    if (is.null(fixedGaps)) {
      # Then we start with a complete graph
      # (i.e,. all nodes connected to each other)
      adjM <- matrix(TRUE, nrow = p, ncol = p)
    }  else if (!identical(dim(fixedGaps), c(p, p))) {
      stop("Dimensions of the dataset and fixedGaps do not agree.")
    } else if (!identical(fixedGaps, t(fixedGaps)) ) {
      stop("fixedGaps must be symmetric")
    } else {
      adjM <- !fixedGaps
    }

    colnames(adjM) <- rownames(adjM) <- labels
    diag(adjM) <- FALSE

    if (any(is.null(fixedEdges))) { ## MM: could be sparse
      fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
    } else if (!identical(dim(fixedEdges), c(p, p))) {
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    } else if (!identical(fixedEdges, t(fixedEdges)) ) {
      stop("fixedEdges must be symmetric")
    }

    ## Creating sepsets
    ####################

    seq_p <- seq_len(p) # creates a sequence of integers from 1 to p
    sepset <- lapply(seq_p, function(.) vector("list",p)) # a list of lists [p x p]

    ##### Initializing variables
    #############################

    mec_score_df <- ord_pag_list_score_df <- diff_citestResults <- c()
    citestResults <- c()

    pag_List <- cur_ord_pag_list <- list()

    if (verbose) {
      cat("Starting with a complete PAG...\n")
    }
    mec_out <- getMEC(adjM * 1, ag.type="pag", scored=TRUE, max.ord = NA,
                      citestResults, indepTest, suffStat,
                      verbose=verbose > 1, allowNewTests=allowNewTests)
    citestResults <- mec_out$citestResults

    cur_ord_pag_list[[length(pag_List)+1]] <-
      list(amat.pag=adjM * 1,
           sepset=sepset,
           sepsetResults=c(),
           doneEdges=fixedEdges | !adjM,
           scores = NULL,
           curord = -1,
           violations=FALSE,
           done=FALSE,
           mec=mec_out$mec,
           ordPAGs=list())

    ord <- 0
    diff_citestResults <- list() # symmetric difference of citests for all candidate PAGs
    error_message <- "No error"
    #library(futile.logger)

    # Set up logging to file
    if (!file.exists(log_folder))
      dir.create(log_folder, recursive = TRUE)

    appender_file <- paste0(file.path(log_folder, "dcfci_log"),
                            format(Sys.time(), '%Y%m%d_%H%M%S%OS.3'), ".txt")
    futile.logger::flog.appender(futile.logger::appender.file(appender_file),
      name = "dcfci_log")
    futile.logger::flog.info(paste("Initializing..."), name = "dcfci_log")


    exceeded_list_max <- FALSE
    while (ord <= min(m.max, p-2)) {

      prev_cur_ord_pag_list <- cur_ord_pag_list

      todo_pag_ids <- c()
      if (length(cur_ord_pag_list) > 0) {
        todo_pag_ids <- which(sapply(cur_ord_pag_list, function(x) {
          any(!x[["doneEdges"]])}))
      }

      if (length(todo_pag_ids) == 0 || length(todo_pag_ids) >= list.max) {
        break
      }

      # todo_pag_ids indicate the PAGs with top scores from previous ord
      todo_pag_list <- cur_ord_pag_list[todo_pag_ids]
      cur_ord_pag_list <-  cur_ord_pag_list[-todo_pag_ids]

      if (verbose)
        cat("length of todo_pag_list in ord=: ", ord, ": ", length(todo_pag_list), "\n")

      toProcessPAGLists <- list()
      length_toProcessPAGLists <- 0
      has_errors <- FALSE
      for (i in seq_along(todo_pag_list))  {
        curpag <- todo_pag_list[[i]]
        # updating mec of curpag to the current ord
        if (verbose)
          cat("Updating MEC skel of previous best candidate PAG", i, "/",
              length(todo_pag_list), " to ord", ord, "...\n")

        skel_out <- updateSkelScore(curpag, ord, citestResults, indepTest,
                                    suffStat, verbose > 1, allowNewTests)
        curpag <- skel_out$apag
        citestResults <- skel_out$citestResults


        if (verbose)
          cat("Getting corresponding list of PAGs to process", i, "/",
              length(todo_pag_list), "...\n")

        # this will potentially run in parallel, so the outer loop
        # cannot run in parallel as well.
        outList <- getToProcessPAGList(
          curpag, ord, fixedEdges, alpha, citestResults,
          indepTest, suffStat, verbose > 1, pH0ThreshMin, pH0ThreshMax, allowNewTests)

        if (!is.null(outList)) {
          toProcessPAGLists[[i]] <- outList
          length_toProcessPAGLists <- length_toProcessPAGLists + length(outList)
          if (verbose)
            cat("curord:", ord, "; current length of toProcessPAGList: ", length_toProcessPAGLists, "\n")
          if (length_toProcessPAGLists > list.max) {
            cat("ERROR: The maximum list size of", list.max, "has been exceeded.\n")
            error_message <- paste0("Reached maximum list size at order ", ord,
                                    " -- list size: ", length_toProcessPAGLists, ".")
            exceeded_list_max <- TRUE
            has_errors <- TRUE
            break
          }
        } else {
          has_errors <- TRUE
          toProcessPAGLists <- NULL
          error_message <- paste0("No PAG to be processed at order ", ord, ".")
          break
        }
      }

      toProcessPAGList <- list()
      if (!has_errors) {
        all_citestResults <- citestResults
        for (i in seq_along(toProcessPAGLists)) {
          if (is.null(toProcessPAGLists[[i]]$toProcessPAGList) || length(toProcessPAGLists[[i]]$toProcessPAGList) == 0) {
            # Adding the PAG without applying any minimal conditional independencies,
            # which requires no processing in this ord
            #if (!all(toProcessPAGLists[[i]]$curpag$doneEdges)) {
              cur_ord_pag_list[[length(cur_ord_pag_list)+1]] <- toProcessPAGLists[[i]]$curpag
            #}
          } else {
            toProcessPAGList <- c(toProcessPAGList,
                                  toProcessPAGLists[[i]]$toProcessPAGList)
          }
          all_citestResults <- rbind(all_citestResults, toProcessPAGLists[[i]]$citestResults)
        }
        # tests from previous order plus those to find separators at the current order.
        citestResults <- all_citestResults[!duplicated(all_citestResults[,1:4]),]

        toProcessPAGList <- getUniquePAGList(toProcessPAGList, key="sepset")$pag_List

        if (verbose)
          cat("curord:", ord, "; length of toProcessPAGList: ", length(toProcessPAGList), "\n")
        if (length(toProcessPAGList) > list.max) {
          cat("ERROR: The maximum list size of", list.max, "has been exceeded.\n")
          error_message <- paste0("Reached maximum list size at order ", ord,
                                  " -- list size: ", length(toProcessPAGList), ".")
          exceeded_list_max <- TRUE
          has_errors <- TRUE
        }
      }

      if (has_errors) {
        toProcessPAGList <- NULL
        cur_ord_pag_list <- prev_cur_ord_pag_list
        break
      }

      #lapply(toProcessPAGList, function(x) {formatSepset(x$sepset)})

      if (length(toProcessPAGList) > 0) {
        if (run_parallel) {
          # Process each PAG by applying the rules of the FCI, checking for
          # violations, and computing the scored MEC
          out_processed <- future_lapply(1:length(toProcessPAGList),
                                         processPotPAG, toProcessPAGList,
                                         ord, citestResults, indepTest, suffStat,
                                         verbose=verbose > 1,
                                         allowNewTests=allowNewTests,
                                         appender_file=appender_file, future.seed=TRUE)

          # Updating used citestResults and newProcessedPAGList
          all_citestResults <- citestResults
          newProcessedPAGList <- list()
          for (i in seq_along(out_processed)) {
            newProcessedPAGList[[i]] <- out_processed[[i]]$potPAG
            all_citestResults <- rbind(all_citestResults, out_processed[[i]]$citestResults)
          }
        } else {
          all_citestResults <- citestResults
          newProcessedPAGList <- list()
          for (i in seq_along(toProcessPAGList)) {
            if (verbose)
              cat("Processing PAG", i, "/", length(toProcessPAGList), "\n")
            out_processed <- processPotPAG(i, toProcessPAGList, ord,
                                           citestResults, indepTest,
                                           suffStat, verbose=verbose > 1,
                                           allowNewTests=allowNewTests,
                                           appender_file = appender_file)
            all_citestResults <- rbind(all_citestResults, out_processed$citestResults)
            newProcessedPAGList[[i]] <- out_processed$potPAG
          }
        }
        citestResults <- all_citestResults[!duplicated(all_citestResults[,1:4]),]


        #lapply(newProcessedPAGList, function(x) {formatSepset(x$sepset)})

        cur_ord_pag_list <- c(cur_ord_pag_list, newProcessedPAGList)


        valid_pags_ids <- which(!sapply(cur_ord_pag_list, function(x) {
          x$violations }))

        if (length(valid_pags_ids) > 0) {
          symm_out <- getSymmDiffCITestResults(cur_ord_pag_list, ord,
                                             verbose=verbose >1)
          diff_citestResults[[as.character(ord)]] <- symm_out$symmDiffCITests
          cur_ord_pag_list <- symm_out$cur_ord_pag_list
        } else {
          error_message <- paste0("No valid PAGs at order ", ord, ".")
          toProcessPAGList <- NULL
          diff_citestResults[[as.character(ord)]] <- data.frame()
          cur_ord_pag_list <- prev_cur_ord_pag_list
          break
        }
      } else {
        diff_citestResults[[as.character(ord)]] <- data.frame()
      }

      #sapply(cur_ord_pag_list, function(x) {x$violations})


      if (nrow(diff_citestResults[[as.character(ord)]]) == 0) {
        for (i in seq_along(cur_ord_pag_list)) {
          # since all pags in the current order do not differ in their mec citests,
          # their union score is simply the mec score...
          mec_score_tuple <- c(cur_ord_pag_list[[i]]$mec$mec_score[1],
                               1 - cur_ord_pag_list[[i]]$mec$mec_mse, # so that the higher, the better
                               cur_ord_pag_list[[i]]$mec$mec_score[2])
          if (!cur_ord_pag_list[[i]]$violations) {
            cur_ord_pag_list[[i]]$scores <- list(ord_symm_diff_score = c(1,1,1),
                                                 ord_union_mec_score = mec_score_tuple,
                                                 mec_score = mec_score_tuple)#cur_ord_pag_list[[i]]$mec$mec_score,
                                                 #mec_mse=cur_ord_pag_list[[i]]$mec$mec_mse,
                                                 #ord_symm_diff_mse=0)
          } else {
            cur_ord_pag_list[[i]]$scores <- list(ord_symm_diff_score = c(0,0,0),
                                                 ord_union_mec_score = c(0,0,0),
                                                 mec_score = c(0,0,0))#,
                                                 #mec_mse=1,
                                                 #ord_symm_diff_mse=1)
          }
          cur_ord_pag_list[[i]]$curord <- ord
          cur_ord_pag_list[[i]]$ordPAGs[[as.character(ord)]] <- cur_ord_pag_list[[i]]
        }
      }

      #sapply(cur_ord_pag_list, function(x) {x$curord})


      ord_pag_list_score_df <- rankPAGList(cur_ord_pag_list, max_ord = ord,
                                            score_type = "pag_list",
                                            order_by_mse = order_by_mse)

      # mec_score is the MEC-targetted pag score
      mec_score_df <- rankPAGList(cur_ord_pag_list, max_ord = ord,
                                  score_type = "mec",
                                  order_by_mse = order_by_mse)

      #lapply(pag_List, function(x) {(x$mec)})
      #lapply(pag_List, function(x) {(x$scores)})

      top_dc_pag_ids <- getTopPagIds(ord_pag_list_score_df,
                                     mec_score_df =mec_score_df,
                                     ord, sel_top, prob_sel_top)

      # subset(mec_score_df, pag_list_id %in% top_dc_pag_ids)
      # subset(ord_pag_list_score_df, pag_list_id %in% top_dc_pag_ids)

      # Adding to the pag_List all processed PAG that have some violation or
      # are not the best of the current order.
      pag_List <- c(pag_List, cur_ord_pag_list[-top_dc_pag_ids])

      #lapply(cur_ord_pag_list, function(x) {formatSepset(x$sepset)})
      #sapply(pag_List, function(x) {x$curord})

      # The other ones proceed to the next order.
      cur_ord_pag_list <- cur_ord_pag_list[top_dc_pag_ids]

      #lapply(cur_ord_pag_list, function(x) {formatSepset(x$sepset)})
      #sapply(cur_ord_pag_list, function(x) {x$curord})

      ord = ord + 1
      if (verbose)
        cat("ord=", ord, "\n")
    }

    pag_List <- c(pag_List, cur_ord_pag_list)
    #lapply(pag_List, function(x) {formatSepset(x$sepset)})
    #sapply(pag_List, function(x) {x$curord})

    top_dcPAGs <- top_scoresDF <- NULL
    if (length(pag_List) > 1) {
      ord_pag_list_score_df <- rankPAGList(pag_List, max_ord = ord-1,
                                           score_type = "pag_list")

      mec_score_df <- rankPAGList(pag_List, max_ord =  ord-1,
                                  score_type = "mec")


      # Putting the score tables again in the order of the pag ids
      ord_pag_list_score_df <- ord_pag_list_score_df[
        order(ord_pag_list_score_df$pag_list_id), ]

      # PAGs are now sorted according to their MEC-Targeted Score (mec_score_df)
      pag_List <- pag_List[mec_score_df$pag_list_id]

      # First tables are sorted according to the order in mec_score_df
      ord_pag_list_score_df <- ord_pag_list_score_df[mec_score_df$pag_list_id, ]


      # Then the order is set sequentially across all score tables
      mec_score_df$pag_list_id <- 1:length(pag_List)
      ord_pag_list_score_df$pag_list_id <- 1:length(pag_List)

      #selects those with top mec_score upper bound at the last iteration
      if (prob_sel_top) {
        top_ids <- which(mec_score_df$prob_index <= prob_sel_top &
                           mec_score_df$violations == FALSE &
                           mec_score_df$duplicated == FALSE)
      } else {
        top_ids <- which(mec_score_df$index <= sel_top &
                           mec_score_df$violations == FALSE &
                           mec_score_df$duplicated == FALSE)
      }
      top_dcPAGs <- pag_List[top_ids]
      top_scoresDF <- mec_score_df[top_ids,]
    }
  })

  return(list(top_dcPAGs=top_dcPAGs, # PAGs ranked as 1st
              top_scoresDF=top_scoresDF, #MEC-targeted scores of top_dcPAGs
              allPAGList=pag_List, # all constructed valid / invalid / duplicated PAGs
              mec_score_df=mec_score_df,
              ord_pag_list_score_df=ord_pag_list_score_df,
              diff_citestResults=diff_citestResults,
              elapsed_time=as.numeric(elapsed_time)[3],
              exceeded_list_max=exceeded_list_max,
              order_processed=ord-1,
              error_message = error_message))
}

getProbIndices <- function(min_scores, max_scores) {
  ind <- 1
  cur_grp <- 1
  grp <- rep(NA, length(min_scores))
  while (ind <= length(min_scores)) {
    cur_min <- round(min_scores[ind], 10)
    pot_next_inds <- which(round(max_scores, 10) < cur_min)
    if (length(pot_next_inds) > 0) {
      next_ind <- pot_next_inds[1]
      grp[ind:(next_ind-1)] <- cur_grp
      ind <- next_ind
    } else {
      grp[ind:length(grp)] <- cur_grp
      break
    }
    cur_grp <- cur_grp + 1
  }
  return(grp)
}

# score_type = "pag_list" or "mec"
# "pag_list" uses the symmetric difference and union
# of the mec scores of all pags in the list
# "mec" is a standalone score, but cannot be used for comparisons.
# score order if score_type == "pag_list":
# ord_symm_diff_up, ord_symm_diff_mse, ord_union_up, ord_union_mse
# ord_symm_diff_low, ord_union_low
rankPAGList <- function(pag_list, max_ord,
                        score_type = "pag_list",
                        order_by_mse = FALSE) {

  score_names=c("ord_symm_diff_score", "ord_union_mec_score")
  score_short_names = c("symmdiff", "union")
  if (score_type == "mec") {
      score_names <- c("mec_score")
      score_short_names = c("mec")
  }

  score_len = 3
  score_labels = c("score_lb", "score_1-mse", "score_up")

  all_scores_list <- list()
  for (score_i in 1:length(score_names)) {
    score_name <- score_names[score_i]
    score_short_name <- score_short_names[score_i]

    scores_df <- getPAGListScores(pag_list, score_name, score_len)

    ord_scores_cols_list <- lapply(scores_df$cur_ord_pag_scores, fromJSON)

    all_scores <- data.frame()
    for (i in seq_along(ord_scores_cols_list)) {
      cur_ord_scores <- ord_scores_cols_list[[i]]
      ord_ints <- c()
      for (j in 0:(max_ord)) {
        cur_int <- (cur_ord_scores[[as.character(j)]])
        if (is.null(cur_int) || length(cur_int) == 0) {
          cur_int <- rep(NA, score_len)
        }
        ord_ints <- c(ord_ints, cur_int)
      }
      all_scores <- rbind.data.frame(all_scores, ord_ints)
    }
    all_scores <- as.data.frame(all_scores)

    # last computed score is the product of all
    colnames(all_scores) <- c(
      paste0("ord", rep(0:max_ord, each=score_len),
             rep(paste0("_", score_short_name, "_", score_labels), max_ord + 1)))

    row.names(all_scores) <- NULL

    all_scores_list[[score_name]] <- all_scores
  }

  all_scores <- all_scores_list[[score_names[1]]]
  if (length(score_names) > 1) {
    for (score_i in 2:length(score_names)) {
      score_name <- score_names[score_i]
      all_scores <- cbind.data.frame(all_scores, all_scores_list[[score_name]])
    }
  }

  top_labels <- c("score_up", "score_1-mse")
  if (order_by_mse) {
    top_labels <- c("score_1-mse", "score_up")
  }

  ordered_labels <- c()
  for (cur_ord in seq(max_ord,0,-1)) {
    for (score_short_name in score_short_names) {
      ordered_labels <- c(ordered_labels,
                        paste0("ord", cur_ord,
                               paste0("_", score_short_name, "_",
                                      top_labels)))
    }
  }
  for (cur_ord in seq(max_ord,0,-1)) {
    for (score_short_name in score_short_names) {
      ordered_labels <- c(ordered_labels,
                          paste0("ord", cur_ord,
                                 paste0("_", score_short_name, "_",
                                        c("score_lb"))))
    }
  }

  tmp_scores <- cbind.data.frame(
    all_scores[, ordered_labels, drop = FALSE],
    pag_list_id=1:nrow(all_scores))

  ordered_ids <- dplyr::arrange_all(tmp_scores, dplyr::desc)$pag_list_id



  tmp_scores[, "violations"] <- scores_df$violation
  tmp_scores <- tmp_scores[ordered_ids, ]

  # all_scores are sorted already...
  index_start_ids <- which(!duplicated(
    tmp_scores[,-which(colnames(tmp_scores) == "pag_list_id")]))
  tmp_scores[index_start_ids, "index"] <- 1:length(index_start_ids)
  for (i in 1:nrow(tmp_scores)) {
    if (is.na(tmp_scores[i, "index"])) {
      tmp_scores[i, "index"] <- tmp_scores[i-1, "index"]
    }
  }

  prob_indices <- getProbIndices(
    tmp_scores[, paste0("ord", max_ord, "_", score_short_names[1], "_score_lb")],
    tmp_scores[, paste0("ord", max_ord, "_", score_short_names[1], "_score_up")])
  new_scores_df <- cbind(tmp_scores,
                         "prob_index"= prob_indices)

  unique_pag_ids <- getUniquePAGList(
    pag_list[new_scores_df$pag_list_id], key="amat")$ids
  new_scores_df$duplicated <- TRUE
  new_scores_df[unique_pag_ids, "duplicated"] <- FALSE

  return(new_scores_df)
}

#
# rankPAGListOLD2 <- function(pag_list, max_ord, scores_df=NULL,
#                         score_name="ord_symm_diff_score",
#                         score_len = 3,
#                         score_labels = c("_score_lb", "_score_1-mse", "_score_up")) {
#
#   if (is.null(scores_df)) {
#     scores_df <- getPAGListScores(pag_list, score_name, score_len)
#   }
#
#   ord_scores_cols_list <- lapply(scores_df$cur_ord_pag_scores, fromJSON)
#   all_scores <- data.frame()
#   for (i in seq_along(ord_scores_cols_list)) {
#     cur_ord_scores <- ord_scores_cols_list[[i]]
#     ord_ints <- c()
#     for (j in 0:(max_ord)) {
#       cur_int <- (cur_ord_scores[[as.character(j)]])
#       if (is.null(cur_int) || length(cur_int) == 0) {
#         cur_int <- rep(NA, score_len)
#       }
#       ord_ints <- c(ord_ints, cur_int)
#     }
#     all_scores <- rbind.data.frame(all_scores, ord_ints)
#   }
#   all_scores <- as.data.frame(all_scores)
#
#
#
#
#   # last computed score is the product of all
#   colnames(all_scores) <- c(
#       paste0("ord", rep(0:max_ord, each=score_len),
#            rep(score_labels, max_ord + 1)))
#
#   row.names(all_scores) <- NULL
#
#   if (grepl("mse", score_name, fixed=TRUE)) { # score_name:  ord_symm_diff_mse" or "mec_mse") {
#     tmp_scores <- cbind.data.frame(
#       all_scores[, c(seq(ncol(all_scores), 1, by=-1)), drop = FALSE],
#       pag_list_id=1:nrow(all_scores))
#     ordered_ids <- dplyr::arrange_all(tmp_scores)$pag_list_id
#   } else {
#     ordered_labels <- c()
#     for (cur_ord in seq(max_ord,0,-1)) {
#       ordered_labels <- c(ordered_labels,
#                           paste0("ord", cur_ord, c("_score_up", "_score_1-mse")))
#     }
#     for (cur_ord in seq(max_ord,0,-1)) {
#       ordered_labels <- c(ordered_labels,
#                           paste0("ord", cur_ord, "_score_lb"))
#     }
#     tmp_scores <- cbind.data.frame(
#       all_scores[, ordered_labels, drop = FALSE],
#       pag_list_id=1:nrow(all_scores))
#
#     ordered_ids <- dplyr::arrange_all(tmp_scores, dplyr::desc)$pag_list_id
#   }
#
#
#   tmp_scores[, "violations"] <- scores_df$violation
#   tmp_scores <- tmp_scores[ordered_ids, ]
#
#   # all_scores are sorted already...
#   index_start_ids <- which(!duplicated(tmp_scores[,-which(colnames(tmp_scores) == "pag_list_id")]))
#   tmp_scores[index_start_ids, "index"] <- 1:length(index_start_ids)
#   for (i in 1:nrow(tmp_scores)) {
#     if (is.na(tmp_scores[i, "index"])) {
#       tmp_scores[i, "index"] <- tmp_scores[i-1, "index"]
#     }
#   }
#
#   if (grepl("mse", score_name, fixed=TRUE)) {
#     new_scores_df <- tmp_scores
#   } else {
#     prob_indices <- getProbIndices(
#       tmp_scores[, paste0("ord", max_ord, "_score_lb")],
#       tmp_scores[, paste0("ord", max_ord, "_score_up")])
#     new_scores_df <- cbind(tmp_scores,
#                            "prob_index"= prob_indices)
#   }
#
#   unique_pag_ids <- getUniquePAGList(
#     pag_list[new_scores_df$pag_list_id], key="amat")$ids
#   new_scores_df$duplicated <- TRUE
#   new_scores_df[unique_pag_ids, "duplicated"] <- FALSE
#
#   return(new_scores_df)
# }
#
#
# rankPAGListOLD <- function(pag_list, max_ord, scores_df=NULL,
#                         score_name="ord_symm_diff_score",
#                         score_len = 2,
#                         score_labels = c("_score_lb", "_score_up")) {
#
#   if (is.null(scores_df)) {
#     scores_df <- getPAGListScores(pag_list, score_name, score_len)
#   }
#
#   ord_scores_cols_list <- lapply(scores_df$cur_ord_pag_scores, fromJSON)
#   all_scores <- data.frame()
#   for (i in seq_along(ord_scores_cols_list)) {
#     cur_ord_scores <- ord_scores_cols_list[[i]]
#     ord_ints <- c()
#     for (j in 0:(max_ord)) {
#       cur_int <- (cur_ord_scores[[as.character(j)]])
#       if (is.null(cur_int) || length(cur_int) == 0) {
#         cur_int <- rep(NA, score_len)
#       }
#       ord_ints <- c(ord_ints, cur_int)
#     }
#     all_scores <- rbind.data.frame(all_scores, ord_ints)
#   }
#   all_scores <- as.data.frame(all_scores)
#
#
#   # last computed score is the product of all
#   colnames(all_scores) <- c(
#     paste0("ord", rep(0:max_ord, each=score_len),
#            rep(score_labels, max_ord + 1)))
#
#   row.names(all_scores) <- NULL
#
#   if (grepl("mse", score_name, fixed=TRUE)) { # score_name:  ord_symm_diff_mse" or "mec_mse") {
#     tmp_scores <- cbind.data.frame(
#       all_scores[, c(seq(ncol(all_scores), 1, by=-1)), drop = FALSE],
#       pag_list_id=1:nrow(all_scores))
#     ordered_ids <- dplyr::arrange_all(tmp_scores)$pag_list_id
#   } else {
#     tmp_scores <- cbind.data.frame(
#       all_scores[, c(seq(ncol(all_scores), 2, by=-2), seq(ncol(all_scores)-1, 1, by=-2)), drop = FALSE],
#       pag_list_id=1:nrow(all_scores)
#     )
#     ordered_ids <- dplyr::arrange_all(tmp_scores, dplyr::desc)$pag_list_id
#   }
#
#
#   tmp_scores[, "violations"] <- scores_df$violation
#   tmp_scores <- tmp_scores[ordered_ids, ]
#
#   # all_scores are sorted already...
#   index_start_ids <- which(!duplicated(tmp_scores[,-which(colnames(tmp_scores) == "pag_list_id")]))
#   tmp_scores[index_start_ids, "index"] <- 1:length(index_start_ids)
#   for (i in 1:nrow(tmp_scores)) {
#     if (is.na(tmp_scores[i, "index"])) {
#       tmp_scores[i, "index"] <- tmp_scores[i-1, "index"]
#     }
#   }
#
#   if (grepl("mse", score_name, fixed=TRUE)) {
#     new_scores_df <- tmp_scores
#   } else {
#     prob_indices <- getProbIndices(
#       tmp_scores[, paste0("ord", max_ord, "_score_lb")],
#       tmp_scores[, paste0("ord", max_ord, "_score_up")])
#     new_scores_df <- cbind(tmp_scores,
#                            "prob_index"= prob_indices)
#   }
#
#   unique_pag_ids <- getUniquePAGList(pag_list, key="amat")$ids
#   new_scores_df$duplicated <- TRUE
#   new_scores_df[unique_pag_ids,"duplicated"] <- FALSE
#
#
#   return(new_scores_df)
# }

getAllPotSepsetsXY <- function(x, y, amat.pag, m.max=Inf, verbose=TRUE) {
  done = FALSE
  ord = 0
  sepset_list <- c()
  while (ord <= m.max && !done) {
    out <- getOrdPotSepsetsXY(x, y, ord, amat.pag, verbose=verbose)
    sepset_list <- c(sepset_list, out$potsepsets)
    done <- out$done
    ord = ord + 1
  }
  return(sepset_list)
}


getUniquePAGList <- function(pag_List, key="sepset") {
  unique_pag_list <- list()
  unique_ids <- c()
  if (key == "sepset") {
    sepsetResultsList <- lapply(pag_List, function(x) {
      ret <- x$sepsetResults
      rownames(ret) <- NULL
      ret
    })
    unique_ids <- which(!duplicated(sepsetResultsList))
  } else { #key = "amat"
    amatList <- lapply(pag_List, function(x) {
      ret <- x$amat.pag
      rownames(ret) <- NULL
      ret
    })
    unique_ids <- which(!duplicated(amatList))
  }
  unique_pag_list <- pag_List[unique_ids]

  return(list(pag_List=unique_pag_list, ids=unique_ids))
}

# sepset contains only one sepset per pair
# sepsets contains a list of sepsets per pair
# the function checks whether each pair's sepset is in the
# respective list of sepsets in the object sepsets.
sepsetInSepsets <- function(sepset, sepsets) {
  p <- length(sepset)
  for (x in 1:p) {
    for (y in 1:p) {
      matched = FALSE
      if (is.null(sepset[[x]][[y]])) {
        matched <- is.null(sepsets[[x]][[y]])
      } else {
        cur_sepsets <- sepsets[[x]][[y]]
        if (!is.list(cur_sepsets)) {
          cur_sepsets <- list(cur_sepsets)
        }
        for (cur_sepset in cur_sepsets) {
          if (isTRUE(all.equal(cur_sepset, sepset[[x]][[y]]))) {
            matched = TRUE
            break
          }
        }
      }
      if (!matched) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#' @export getSepsetIndex
getSepsetIndex <- function(pag_List, sepset) {
  sepset_List <- lapply(pag_List, function(x) {x$sepset})
  index <- NA
  for (i in seq_along(sepset_List)) {
    if (sepsetInSepsets(sepset_List[[i]], sepset)) {
      index <- i
      break
    }
  }
  return(index)
}

