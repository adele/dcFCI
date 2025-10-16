#' @param amat.ag: the pcalg adjancency matrix of an ancestral graph
#' @param ag.type: the ancestral graph's type, either 'pag' or 'mag'
#' @importFrom pcalg pcalg2dagitty
#' @importFrom FCI.Utils MAGtoMEC updateMECTripletIsolators getMAGImpliedSepset
#' @export getMEC
getMEC <- function(amat.ag, ag.type="pag", scored=FALSE, max.ord=NA,
                   citestResults=NULL, indepTest=NULL, suffStat=NULL,
                   verbose=FALSE, allowNewTests=TRUE) {

  if (!isValidPAG(amat.ag)) {
    return(NULL)
  }

  if (ag.type == "mag") {
    amat.mag <- amat.ag
  } else if (ag.type == "pag") {
    mag_out <- getMAG(amat.ag)
    magg <- mag_out$magg
    amat.mag <- mag_out$amat.mag
    amat.mag <- amat.mag[colnames(amat.ag), colnames(amat.ag)]
  } else {
    stop("ag.type must be \'pag\' or \'mag\'")
  }

  labels <- colnames(amat.ag)
  mec <- FCI.Utils::MAGtoMEC(amat.mag, verbose)
  mec <- FCI.Utils::updateMECTripletIsolators(mec, amat.ag, verbose)

  if (scored) {
    magg <- pcalg::pcalg2dagitty(amat.mag, labels, type = "mag")
    sepset <- FCI.Utils::getMAGImpliedSepset(magg, labels)

    if (is.na(max.ord)) {
      max.ord <- getSepsetMaxOrd(sepset)
      if (is.na(max.ord)) {
        # then it is a complete graph, with no user-defined max.ord
        # assuming that the discovery has not started yet
        max.ord <- 0
      }
    }
    scores_out <- scoreMEC(mec, sepset, max.ord, citestResults, indepTest,
                           suffStat, verbose, allowNewTests)
    citestResults <- scores_out$citestResults
    mec <- scores_out$scored_mec
  }
  return(list(mec=mec,
              citestResults=citestResults))
}
