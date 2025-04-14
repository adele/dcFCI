sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}


#' @importFrom FCI.Utils getMConnPaths isCollider isDefiniteNonCollider
#' @importFrom pcalg searchAM
#' @export getImpliedConditionalSepset
getImpliedConditionalSepset <- function(pagAdjM, xname, yname, snames, definite=TRUE,
                                        ignore_path_list = list(),
                                        ignore_sepvar_names = list(),
                                        verbose = FALSE) {
  Z <- snames
  updated = TRUE
  while (updated) {
    defconnpaths <- getMConnPaths(pagAdjM, xname, yname, Z, definite = definite)
    ndcpaths <- length(defconnpaths)
    rmids <- c()
    if (ndcpaths > 0) {
      for (i in 1:ndcpaths) {
        dcpath <- defconnpaths[[i]]
        if (sapply(ignore_path_list, function(x) {
          length(x) == length(dcpath) && all(x %in% dcpath)} )) {
          rmids <- c(rmids, i)
        }
      }
    }
    pathList <- defconnpaths
    if (!is.null(rmids)) {
      pathList <- defconnpaths[-rmids]
    }
    i = 1
    updated = FALSE
    for (curpath in pathList) {
      # there are still paths that are connecting and of definite status
      len_curpath <- length(curpath)

      if (len_curpath > 2) {
        for (k in 1:(len_curpath-2)) {
          vi_name <- curpath[k]
          vm_name <- curpath[k+1]
          vj_name <- curpath[k+2]

          triplet <- c(vi_name, vm_name, vj_name)
          vnames <- colnames(pagAdjM)
          if (isCollider(pagAdjM, triplet)) {
            vm_id <- which(vnames == triplet[2])
            devm <- vnames[pcalg::searchAM(pagAdjM, vm_id, type="de")]
            if (length(which(devm %in% Z)) == 0) {
              # no definite descendant of this collider is in Z, so the path is
              # blocked and we don't need to checked it again util Z is changed...
              #pathList <- pathList[-i]
              break
            } else {
              next
            }
          } else if (isDefiniteNonCollider(pagAdjM, triplet)) {
            if ((vm_name %in% ignore_sepvar_names)) {
              next
            } else {
              if (length(which(Z == vm_name)) == 0) {
                # we add the non-collider to Z
                Z <- c(Z, vm_name)
                updated = TRUE # since we added a new variable to Z, all paths have to be rechecked.
              }
              # since the path is definitely blocked, we go to the next path in pathList
              break
            }
          }
        }
      }
    }
  }
  return(list(sepset=Z, connpaths = pathList))
}

#' @export getCommonPAG
getCommonPAG <- function(amatPagList) {
  if (is.null(amatPagList) || length(amatPagList) < 1) {
    return(NULL)
  }

  common_pag <- amatPagList[[1]]
  for (i in 1:length(amatPagList)) {
    common_pag[which(common_pag - amatPagList[[i]] != 0, arr.ind = T)] <- 1
  }
  return(common_pag)
}

#' @export getSepsetResults
getSepsetResults <- function(citestResults, sepset) {
  sepsetResults <- c()
  p <- length(sepset)
  for (i in 1:p) {
    for (j in i:p) {
      cur_sepsets <- sepset[[i]][[j]]
      if (!is.null(cur_sepsets)) {
        if (!is.list(cur_sepsets)) {
          cur_sepsets <- list(cur_sepsets)
        }
        for (cursepset in cur_sepsets) {
          curSStr <- getSepString(cursepset)
          curResults <- subset(citestResults, X==i & Y==j & S==curSStr)
          if (length(curResults) > 0) {
            sepsetResults <- rbind(sepsetResults,curResults)
          } else {
            warning("no citestResult for X=", x, "Y=", y, "S=", curSStr)
          }
        }
      }
    }
  }
  return(sepsetResults)
}

#' @export getEssentialCITestResults
getEssentialCITestResults <- function(citestResults, sepset, max.ord=NULL) {
  if (is.null(max.ord)) {
    max.ord <- length(sepset) - 2
  }

  sepsetResults <- c()
  p <- length(sepset)
  for (i in 1:p) {
    for (j in i:p) {
      cur_sepsets <- sepset[[i]][[j]]
      if (!is.null(cur_sepsets)) {
        if (!is.list(cur_sepsets)) {
          cur_sepsets <- list(cur_sepsets)
        }
        for (cursepset in cur_sepsets) {
          cur_subsets <- getSubsets(cursepset, only_proper = FALSE)
          for (cur_subset in cur_subsets) {
            curSStr <- getSepString(cur_subset)
            curResults <- subset(citestResults, X==i & Y==j & S==curSStr)
            if (length(curResults) > 0) {
              sepsetResults <- rbind(sepsetResults, curResults)
            } else {
              warning("no citestResult for X=", x, "Y=", y, "S=", curSStr)
            }
          }
        }
      } else {
        curResults <- subset(citestResults, X==i & Y==j & ord <= max.ord)
        sepsetResults <- rbind(sepsetResults, curResults)
      }
    }
  }
  return(sepsetResults)
}

#' @export getRPAG
getRPAG <- function(apag, r, sepset=NULL, verbose=FALSE) {
  cur_r_sepset <- if (is.null(sepset)) getPAGImpliedSepset(apag) else sepset
  cur_r_skel <- apag > 0
  p <- length(cur_r_sepset)
  for (i in 1:p) {
    for (j in i:p) {
      cur_sepsets <- cur_r_sepset[[i]][[j]]
      if (!is.null(cur_sepsets)) {
        if (is.list(cur_sepsets) & length(cur_sepsets) > 0) {
          if (length(cur_r_sepset[[i]][[j]][[1]]) == 0) {
            cur_r_sepset[[i]][[j]] <- integer(0)
            cur_r_sepset[[j]][[i]] <- integer(0)
          } else {
            cur_r_sepset[[i]][[j]] <- cur_r_sepset[[i]][[j]][[1]] # using pcalg sepset format
            cur_r_sepset[[j]][[i]] <- cur_r_sepset[[i]][[j]] # using pcalg sepset format
          }
        }
        cur_sepsets <- cur_r_sepset[[i]][[j]]
        if (length(cur_sepsets) > r) {
          cur_r_sepset[[i]][j] <- list(NULL)
          cur_r_sepset[[j]][i] <- list(NULL)
          cur_r_skel[i,j] = cur_r_skel[j,i] <- 1
        }
      }
    }
  }
  #print(formatSepset(cur_r_sepset))

  rules <- rep(TRUE, 10)
  cur_r_amat.pag <- pcalg::udag2pag(pag = cur_r_skel, cur_r_sepset, rules = rules,
                                  unfVect = c(), verbose = verbose > 1, orientCollider = TRUE)
  #renderAG(cur_r_amat.pag)
  return(list(r.amat.pag = cur_r_amat.pag, r.sepset = cur_r_sepset))
}
