getDCFCIMetrics <- function(dcfci_out, dat, citestResults, true.amat.pag) {
  top_mec_score_up <- dcfci_out$mec_score_df[1,1]
  top_mec_score_1mse <- dcfci_out$mec_score_df[1,2]

  index1_pags <- dcfci_out$allPAGList[
    which(dcfci_out$mec_score_df[,1] >= top_mec_score_up &
          dcfci_out$mec_score_df[,2] >= top_mec_score_1mse &
          dcfci_out$mec_score_df$violations == FALSE &
          dcfci_out$mec_score_df$duplicated == FALSE)]

  probindex1_pags <- dcfci_out$allPAGList[
    which(dcfci_out$mec_score_df$prob_index == 1 &
            dcfci_out$mec_score_df$duplicated == FALSE)]

  ntop <- length(index1_pags)
  nprobtop <- length(probindex1_pags)


  truePAGListInd <- getSepsetIndex(dcfci_out$allPAGList, true.sepset)
  truePAGInd <- truePAGProbInd <- NA
  if (!is.na(truePAGListInd)) {
    truePAGInd <- dcfci_out$mec_score_df[truePAGListInd, "index"]
    truePAGProbInd <- dcfci_out$mec_score_df[truePAGListInd, "prob_index"]
  }

  cur_dcfci_metrics <- data.frame()
  for (dcfci_top_pag in index1_pags) {
    est.amat.pag <- dcfci_top_pag$amat.pag
    #renderAG(est.amat.pag)
    est.sepset <- dcfci_top_pag$sepset
    #print(formatSepset(est.sepset))

    dcfci_metrics_partial <-  getMetrics(true.amat.pag, est.amat.pag, est.sepset,
                                         dat=dat, conservative = FALSE)

    if (is.null(dcfci_top_pag$mec$mec_mse)) {
      mec_probs <- c(subset(dcfci_top_pag$mec$all_citests, type == "indep")$pH0,
                     subset(dcfci_top_pag$mec$all_citests, type == "dep")$pH1)
      dcfci_top_pag$mec$mec_mse <- getNormalizedSquaredL2DistanceFromCertainty(mec_probs)
    }

    trivial_score <- getStraightforwardPAGScore(est.amat.pag, suffStat$citestResults)
    cur_dcfci_metrics <- rbind.data.frame(cur_dcfci_metrics,
                                          cbind(dcfci_metrics_partial,
                                                mec_score = t(dcfci_top_pag$mec$mec_score), # not comparable
                                                mec_score_mse = dcfci_top_pag$mec$mec_mse,
                                                trivial_score = t(trivial_score$score), # comparable
                                                trivial_score_mse = trivial_score$mse
                                          ))
  }

  cur_dcfci_metrics_min <- cur_dcfci_metrics
  if (dim(cur_dcfci_metrics)[1] > 1) {
    cur_dcfci_metrics_min <- as.data.frame(t(sapply(cur_dcfci_metrics, min)))
    cur_dcfci_metrics <- as.data.frame(t(colMeans(cur_dcfci_metrics)))
  }


  time_taken <- if (length(dcfci_out$elapsed_time) == 1) {
      as.numeric(dcfci_out$elapsed_time)
    } else {
      as.numeric(dcfci_out$elapsed_time)[3]
    }

  cur_dcfci_metrics <- cbind.data.frame(
    cur_dcfci_metrics,
        truePAGInd=truePAGInd,
        truePAGProbInd = truePAGProbInd,
        max_ord =  dcfci_out$order_processed,
        max_reached = dcfci_out$exceeded_list_max,
        time_taken = time_taken,
        npags = length(dcfci_out$allPAGList),
        ntop = ntop,
        nprobtop=nprobtop)

  cur_dcfci_metrics_min <- cbind.data.frame(
    cur_dcfci_metrics_min,
    truePAGInd=truePAGInd,
    truePAGProbInd = truePAGProbInd,
    max_ord =  dcfci_out$order_processed,
    max_reached = dcfci_out$exceeded_list_max,
    time_taken = time_taken,
    npags = length(dcfci_out$allPAGList),
    ntop = ntop,
    nprobtop=nprobtop)


  return(list(dcfci_metrics=cur_dcfci_metrics,
              dcfci_metrics_min=cur_dcfci_metrics_min))
}
