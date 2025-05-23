getOrdPotSepsets <- function(pagAdjM, ord, alpha, amec,
                             fixedEdges, doneEdges, citestResults,
                             indepTest, suffStat, NAdelete, verbose,
                             allowNewTests=TRUE) {
  UnFaithful = S = NULL
  bestSepResults <- c()

  adjM <- pagAdjM > 0
  # indices of the pairs of connected variables
  ind <- as.data.frame(which(lower.tri(adjM) & adjM, arr.ind = TRUE))
  rownames(ind) <- NULL

  # iterates over all edges with both edgemarks not necessarily faithful and then
  # edges including faithful edgemarks
  remEdge_i <- 1
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
                              suffStat, NAdelete, verbose, allowNewTests)
          pvalue <- citest_out$pvalue
          pH0 <- citest_out$pH0
          pH1 <- citest_out$pH1
          citestResults <- citest_out$citestResults

          if (pvalue > alpha || pH0 > 0.5) { # independent
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

  # We now sort using pvalue
  if (!is.null(bestSepResults) && nrow(bestSepResults) > 0) {
    bestSepResults <- bestSepResults[order(bestSepResults$pvalue, decreasing = T),]
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
getPAGListScores <- function(pag_List, ord_global_score="global_score") {
  violations <- sapply(pag_List, function (x) {x$violations})

  ord_global_scores <- data.frame()
  for (pag in pag_List) {
    cur_ord_global_scores <- as.character(toJSON(
      lapply(pag$ordPAGs, function(x) {
        if (!is.null(x$scores)) {
          x$scores[[ord_global_score]]
        } else {
          rep(0, 2)
        }
      }), digits=10))
    ord_global_scores <- rbind.data.frame(ord_global_scores, data.frame(cur_ord_global_scores))
  }

  scores_df = cbind.data.frame(violation=violations, ord_global_scores)
  colnames(scores_df) <- c("violation", "cur_ord_global_scores")

  return(scores_df)
}


updateSkelScore <- function(apag, ord, citestResults, indepTest,
                            suffStat, NAdelete, verbose,
                            allowNewTests) {
  out_skel <- scoreSkel(apag$sepset, ord,
                        citestResults, indepTest,
                        suffStat, NAdelete, verbose > 1,
                        allowNewTests=allowNewTests)
  apag$mec$skel_citests <- out_skel$skel_citests
  apag$mec$skel_score <- out_skel$skel_score
  apag$mec$all_citests <- unique(rbind(
    apag$mec$all_citests, apag$mec$skel_citests))

  apag$mec$mec_score <- getProbConjunction(
    c(subset(apag$mec$all_citests, type == "indep")$pH0,
      subset(apag$mec$all_citests, type == "dep")$pH1))

  citestResults <- out_skel$citestResults

  return(list(apag=apag, citestResults=citestResults))
}


processPotPAG <- function(potPAG, ord, citestResults, indepTest, suffStat, NAdelete,
                          verbose=FALSE, allowNewTests=TRUE) {
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
                      indepTest, suffStat, NAdelete,
                      verbose=verbose > 1, allowNewTests=allowNewTests)
    potPAG$mec <- out_mec$mec
    citestResults <- out_mec$citestResults
  } else {
    potPAG$mec <- NULL
  }

  potPAG$done <- TRUE

  return(list(potPAG=potPAG, citestResults=citestResults))
}


#' @importFrom rje powerSet
#' @importFrom jsonlite fromJSON
#' @importFrom future.apply future_lapply
#' @importFrom pcalg udag2pag
#' @export dcFCI
dcFCI <- function(suffStat, indepTest, labels, alpha=0.05,
                  NAdelete = TRUE, m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                  verbose = TRUE, sel_top = 1, prob_sel_top = FALSE,
                  run_parallel = TRUE, allowNewTests=TRUE,
                  list.max = 500) {

  score_type = "cur_diffs"
  ord_global_score = "global_score"

  if (is.null(prob_sel_top)) {
    prob_sel_top = 1
  }

  if (is.null(sel_top)) {
    sel_top = 1
  }

  getToProcessPAGList <- function(curpag, ord, fixedEdges, alpha, citestResults,
                                  indepTest, suffStat, NAdelete, verbose,
                                  allowNewTests=TRUE) {
    pagAdjM <-  curpag$amat.pag
    p <- ncol(pagAdjM)
    doneEdges <- curpag$doneEdges
    sepset <- curpag$sepset
    amec <- curpag$amec
    sepsetResultsOut <- getOrdPotSepsets(pagAdjM, ord, alpha, amec,
                                         fixedEdges, doneEdges, citestResults,
                                         indepTest, suffStat, NAdelete, verbose,
                                         allowNewTests=allowNewTests)

    citestResults <- sepsetResultsOut$citestResults
    sepsetResults <- sepsetResultsOut$sepsetResults

    # updating the input pag (i.e, curpag) in the output list.
    curpag$doneEdges <- sepsetResultsOut$doneEdges
    curpag$done <- TRUE # its sepset won't change

    # this creates a pag corresponding to each possible combination of the sepsets
    toProcessPAGList <- list()
    if (!is.null(sepsetResults) && nrow(sepsetResults) > 0) {
      unique_sep_pairs <- unique(sepsetResults[,c(2:3)])
      listS <- list()
      for (pair_row in 1:nrow(unique_sep_pairs)) {
        listS[[pair_row]] <- subset(sepsetResults, X==unique_sep_pairs[pair_row,1] &
                                  Y == unique_sep_pairs[pair_row,2], select = "S", drop=TRUE)
      }


      allCombSepsetPairs <- powerSet(1:nrow(unique_sep_pairs))[-1]  # removes the empty set

      for (k in 1:length(allCombSepsetPairs)) {
        curSepPairs <- allCombSepsetPairs[[k]]
        if (verbose) {
          print(paste0("Generating toProcessPAGList for curSepPairs number: ", k, "/", length(allCombSepsetPairs)))
        }
        # pairs as rows of the unique_sep_pairs
        combinationsS <- expand.grid(listS[curSepPairs], stringsAsFactors = FALSE)

        curToProcessPAGList <-
          foreach (combSrow = 1:nrow(combinationsS),
                   .options.future = list(seed = TRUE))  %dofuture% {
            potAdjM <- curpag$amat.pag > 0
            potSepset <- sepset
            doneEdges <- curpag$doneEdges
            curSepsetResults <- c()

            # iterating over the current pair rows of the unique_sep_pairs
            for (i in 1:length(curSepPairs)) {
              pair_rowid <- curSepPairs[i]
              x = unique_sep_pairs[pair_rowid,1]
              y = unique_sep_pairs[pair_rowid,2]
              Sstr = combinationsS[combSrow, i]

              S <- FCI.Utils::getSepVector(Sstr)
              curSepsetResults <- rbind(curSepsetResults,
                                        subset(sepsetResults, X == x & Y == y & S == Sstr))

              potAdjM[x, y] <- potAdjM[y, x] <- 0
              potSepset[[x]][[y]] <- potSepset[[y]][[x]] <- S
              doneEdges[x, y] <- doneEdges[y, x] <- TRUE
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
        toProcessPAGList <- c(toProcessPAGList, curToProcessPAGList)
      }
    }
    return(list(curpag=curpag,
                toProcessPAGList=toProcessPAGList,
                citestResults=citestResults))
  }

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

  citestResults <- c()

  pag_List <- cur_ord_pag_list <- list()

  if (verbose) {
    cat("Starting with a complete PAG...\n")
  }
  mec_out <- getMEC(adjM * 1, ag.type="pag", scored=TRUE, max.ord = NA,
                    citestResults, indepTest, suffStat, NAdelete,
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

  #system.time({
  while (ord <= min(m.max, p-2)) {
    prev_cur_ord_pag_list <- cur_ord_pag_list

    todo_pag_ids <- which(sapply(cur_ord_pag_list, function(x) {
        any(!x[["doneEdges"]])}))

    if (length(todo_pag_ids) == 0 || length(todo_pag_ids) >= list.max) {
      break
    }

    # todo_pag_ids indicate the PAGs with top scores from previous ord
    todo_pag_list <- cur_ord_pag_list[todo_pag_ids]
    cur_ord_pag_list <-  cur_ord_pag_list[-todo_pag_ids]

    if (verbose)
      cat("length of todo_pag_list in ord=: ", ord, ": ", length(todo_pag_list), "\n")

    if (run_parallel) {
      # For each of the current ord's PAG, we generate
      # all possible combinations of new sepsets
      toProcessPAGLists <- foreach (i = 1:length(todo_pag_list),
           .options.future = list(seed = TRUE)) %dofuture% {
             curpag <- todo_pag_list[[i]]
             # updating mec of curpag to the current ord
             if (verbose)
               cat("Updating MEC skel of previous best candidate PAG", i, "to ord", ord, "...\n")

             skel_out <- updateSkelScore(curpag, ord, citestResults, indepTest,
                             suffStat, NAdelete, verbose > 1,
                             allowNewTests)
             curpag <- skel_out$apag

             getToProcessPAGList(
               curpag, ord, fixedEdges, alpha, skel_out$citestResults,
               indepTest, suffStat, NAdelete, verbose > 1, allowNewTests)
           }
    } else {
      toProcessPAGLists <- list()
      for (i in 1:length(todo_pag_list))  {
        curpag <- todo_pag_list[[i]]
        # updating mec of curpag to the current ord
        if (verbose)
          cat("Updating MEC skel of previous best candidate PAG", i, "to ord", ord, "...\n")

        skel_out <- updateSkelScore(curpag, ord, citestResults, indepTest,
                                    suffStat, NAdelete, verbose > 1,
                                    allowNewTests)
        curpag <- skel_out$apag
        citestResults <- skel_out$citestResults

        toProcessPAGLists[[i]] <- getToProcessPAGList(
               curpag, ord, fixedEdges, alpha, citestResults,
               indepTest, suffStat, NAdelete, verbose > 1, allowNewTests)
       }
    }



    toProcessPAGList <- list()
    all_citestResults <- citestResults
    for (outList in toProcessPAGLists) {
      toProcessPAGList <- c(toProcessPAGList,
                            outList$toProcessPAGList)
      cur_ord_pag_list[[length(cur_ord_pag_list)+1]] <- outList$curpag
      all_citestResults <- rbind(all_citestResults, outList$citestResults)
    }
    # tests from previous order plus those to find separators at the current order.
    citestResults <- all_citestResults[!duplicated(all_citestResults[,1:4]),]

    toProcessPAGList <- getUniquePAGList(toProcessPAGList, key="sepset")$pag_List

    if (verbose)
      cat("curord:", ord, "; length of toProcessPAGList: ", length(toProcessPAGList), "\n")
    if (length(toProcessPAGList) > list.max) {
      cur_ord_pag_list <- prev_cur_ord_pag_list
      break
    }

    if (length(toProcessPAGList) > 0) {
      if (run_parallel) {
        # Process each PAG by applying the rules of the FCI, checking for
        # violations, and computing the scored MEC
        out_processed <- future_lapply(toProcessPAGList, processPotPAG,
                                       ord, citestResults, indepTest, suffStat,
                                       NAdelete, verbose=verbose > 1,
                                       allowNewTests=allowNewTests, future.seed=TRUE)
        all_citestResults <- citestResults
        newProcessedPAGList <- list()
        for (i in 1:length(out_processed)) {
          newProcessedPAGList[[i]] <- out_processed[[i]]$potPAG
          all_citestResults <- rbind(all_citestResults, out_processed[[i]]$citestResults)
        }
      } else {
        all_citestResults <- citestResults
        newProcessedPAGList <- list()
        for (i in 1:length(toProcessPAGList)) {
          if (verbose)
            cat("Processing PAG", i, "/", length(toProcessPAGList), "\n")
          out_processed <- processPotPAG(toProcessPAGList[[i]], ord,
                                         citestResults, indepTest,
                                         suffStat, NAdelete, verbose=verbose > 1,
                                         allowNewTests=allowNewTests)
          all_citestResults <- rbind(all_citestResults, out_processed$citestResults)
          newProcessedPAGList[[i]] <- out_processed$potPAG
        }
      }
      citestResults <- all_citestResults[!duplicated(all_citestResults[,1:4]),]


      cur_ord_pag_list <- c(cur_ord_pag_list, newProcessedPAGList)


      symm_out <- getSymmDiffCITestResults(cur_ord_pag_list, ord,
                                           verbose=verbose >1)

      #lapply(symm_out$cur_ord_pag_list, function(x) {x$scores})

      diff_citestResults[[as.character(ord)]] <- symm_out$symmDiffCITests
      cur_ord_pag_list <- symm_out$cur_ord_pag_list

    } else {
      diff_citestResults[[as.character(ord)]] <- data.frame()
    }

    if (nrow(diff_citestResults[[as.character(ord)]]) == 0) {
      for (i in seq_along(1:length(cur_ord_pag_list))) {
        if (!cur_ord_pag_list[[i]]$violations) {
          cur_ord_pag_list[[i]]$scores <- list(global_score=c(1,1),
                                               mec_score=cur_ord_pag_list[[i]]$mec$mec_score)
        } else {
          cur_ord_pag_list[[i]]$scores <- list(global_score=c(0,0),
                                               mec_score=c(0,0))
        }
        cur_ord_pag_list[[i]]$curord <- ord
        cur_ord_pag_list[[i]]$ordPAGs[[as.character(ord)]] <- cur_ord_pag_list[[i]]
      }
    }

    # score considering the symmdiff of the union
    global_score_df <- rankPAGList(cur_ord_pag_list, max_ord = ord,
                                   score_name = "global_score")
    mec_score_df <- rankPAGList(cur_ord_pag_list, max_ord = ord,
                                   score_name = "mec_score")

    cur_ret <- list(cur_ord_pag_list=cur_ord_pag_list,
                    pag_List=pag_List,
                    citestResults=citestResults,
                    diff_citestResults=diff_citestResults)

    #lapply(pag_List, function(x) {(x$mec)})
    #lapply(pag_List, function(x) {(x$scores)})

    top_dc_pag_ids <- getTopPagIds(global_score_df, mec_score_df, ord, sel_top, prob_sel_top)$top_dc_pag_ids

    # Adding to the pag_List all processed PAG that have some violation or
    # are not the best of the current order.
    pag_List <- c(pag_List, cur_ord_pag_list[-top_dc_pag_ids])

    # The other ones proceed to the next order.
    cur_ord_pag_list <- cur_ord_pag_list[top_dc_pag_ids]

    ord = ord + 1
    if (verbose)
      cat("ord=", ord, "\n")
  }

  pag_List <- c(pag_List, cur_ord_pag_list)


  global_score_df <- rankPAGList(pag_List, max_ord =  ord-1, score_name = "global_score")
  mec_score_df <- rankPAGList(pag_List, max_ord =  ord-1, score_name = "mec_score")
  top_dc_pag_ids_out <- getTopPagIds(global_score_df, mec_score_df, ord-1, sel_top, prob_sel_top)
  if (!(length(top_dc_pag_ids_out$top_global_pags_ids) == length(top_dc_pag_ids_out$top_mec_pags_ids) &&
      all(top_dc_pag_ids_out$top_global_pags_ids %in% top_dc_pag_ids_out$top_mec_pags_ids))) {
    top_dc_pag_ids <- top_dc_pag_ids_out$top_global_pags_ids
    output_top_score <- "global_score"
  } else {
    top_dc_pag_ids <- top_dc_pag_ids_out$top_dc_pags_ids
    output_top_score <- "mec_global_score"
  }

  # Here, we sort the lists and pags according to the global_score
  mec_score_df <- mec_score_df[order(mec_score_df$pag_list_id), ]
  mec_score_df <- mec_score_df[global_score_df$pag_list_id, ]

  pag_List <- pag_List[global_score_df$pag_list_id]
  global_score_df$pag_list_id <- 1:length(pag_List)
  mec_score_df$pag_list_id <- 1:length(pag_List)

  if (prob_sel_top) {
    top_dcPAGs <- pag_List[which(global_score_df$prob_index == 1)]
    top_scoresDF <- subset(global_score_df, prob_index == 1)
  } else {
    top_dcPAGs <- pag_List[which(global_score_df$index == 1)]
    top_scoresDF <- subset(global_score_df, index == 1)
  }

  return(list(top_dcPAGs=top_dcPAGs, # PAGs ranked as 1st
              top_scoresDF=top_scoresDF, #scores of top_dcPAGs
              allPAGList=pag_List, # all constructed valid / invalid / duplicated PAGs
              scoresDF=global_score_df, # scores of all constructed valid / duplicated invalid PAGs
              mec_score_df=mec_score_df,
              citestResults=citestResults,
              diff_citestResults=diff_citestResults))
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

rankPAGList <- function(pag_list, max_ord, scores_df=NULL,
                        score_name="global_score") {

  if (is.null(scores_df)) {
    scores_df <- getPAGListScores(pag_list, score_name)
  }

  ord_scores_cols_list <- lapply(scores_df$cur_ord_global_scores, fromJSON)
  all_scores <- c()
  for (i in 1:length(ord_scores_cols_list)) {
    cur_ord_scores <- ord_scores_cols_list[[i]]
    ord_ints <- c()
    for (j in 0:(max_ord)) {
      cur_int <- cur_ord_scores[[as.character(j)]]
      if (is.null(cur_int)) {
        cur_int <- c(0,0)
      }
      ord_ints <- c(ord_ints, cur_int)
    }
    all_scores <- rbind(all_scores, ord_ints)
  }
  all_scores <- as.data.frame(all_scores)


  # last computed score is the product of all
  colnames(all_scores) <- c(
    paste0("ord", rep(0:max_ord, each=2),
           rep(c("_score_lb", "_score_up"), max_ord + 1)))#,
  #"global_score_lb", "global_score_up")
  row.names(all_scores) <- NULL


  tmp_scores <- cbind(
    all_scores[, c(seq(ncol(all_scores), 2, by=-2), seq(ncol(all_scores)-1, 1, by=-2))],
    pag_list_id=1:nrow(all_scores)
  )

  ordered_ids <- dplyr::arrange_all(tmp_scores, dplyr::desc)$pag_list_id

  tmp_scores[, "violations"] <- scores_df$violation
  tmp_scores <- tmp_scores[ordered_ids, ]

  # all_scores are sorted already...
  index_start_ids <- which(!duplicated(tmp_scores[,-which(colnames(tmp_scores) == "pag_list_id")]))
  tmp_scores[index_start_ids, "index"] <- 1:length(index_start_ids)
  for (i in 1:nrow(tmp_scores)) {
    if (is.na(tmp_scores[i, "index"])) {
      tmp_scores[i, "index"] <- tmp_scores[i-1, "index"]
    }
  }

  new_scores_df <- cbind(tmp_scores,
                         "prob_index"= getProbIndices(
                           tmp_scores[, paste0("ord", max_ord, "_score_lb")],
                           tmp_scores[, paste0("ord", max_ord, "_score_up")]))

  unique_pag_ids <- getUniquePAGList(pag_list, key="amat")$ids
  new_scores_df$duplicated <- TRUE
  new_scores_df[unique_pag_ids,"duplicated"] <- FALSE


  return(new_scores_df)
}


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
  } else {
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

