## Little Library to perform various a posteriori model selection ##

###########################################################################
# compute_AUC -------------------------------------------------------------

#' Compute AUC Index
#'
#' Function to compute the area under ROC curve exploiting trapezoidal rule.
#'
#' @param x Numeric vector of length \code{n}
#' @param y Numeric vector of length \code{n}
#' @return A numeric value which is the area under the ROC curve.

#' @export
compute_AUC = function(x,y) {
  l = length(x)
  if(length(x) != length(y) )
    stop('length of x and y shoud be the same.')
  auc = 0
  for(i in 2:l){
    auc = auc + (x[i] - x[i-1])*(y[i-1] + y[i])/2
  }
  return(auc)
}

###########################################################################

###########################################################################
# misclass_analysis ----------------------------------------------------

#' Bayesian Sensitivity Analysis
#'
#' Given the true graph and a plinks matrix, this utility computes the confusion matrices, plot the ROC curve (if required)
#' and the AUC index, selects the best threshold according to the number of misclassified links
#' and also return the best graph according to the abovementioned analysis.
#' @param PL The posterior probability of inclusion of each link according to \code{\link{FGM_sampler}}.
#' It should be an upper diagonal matrix.
#' @param G_true The true graph.
#' @param tol The sequence of thresholds, used to select a graph truncating PL at that value.
#' @param ROC Boolean parameter. If \code{TRUE}, the plot the ROC curve is showed.
#' @return A list of 5 elements:
#' \itemize{
#'   \item{\strong{all_confusion}, a list with all the confusion matrices, one for each tol value.}
#'   \item{\strong{best_threshold}, the best value of \code{tol} according to this analysis.}
#'   \item{\strong{best_confusion}, the best confusion matrix, the one corresponding to \code{best_threshold}.}
#'   \item{\strong{best_truncated_graph}, the proposed posterior graph according to the analysis.}
#'   \item{\strong{AUC}, the AUC Index corresponding to the best value selected.}
#' }

#' @export
misclass_analysis = function(PL, G_true, tol = seq(0.1, 1, by = 0.05), ROC = FALSE){
  p = dim(PL)[1]
  PL = PL + t(PL)
  diag(PL) <- diag(G_true) <- rep(0,p)
  confusion = list()
  if(any(tol > max(PL)))
    tol <- tol[-which(tol > max(PL))]
  if(is.null(tol))
    stop("No feasible tolerances")
  giusti = rep(0,length(tol))
  sensitivity = rep(0,length(tol))
  uno_meno_specificity = rep(0,length(tol))
  for (i in 1:length(tol)) {
    tolerance <- tol[i]
    Estimated = matrix(0,p,p)
    Estimated[PL >= tolerance] = 1
    confusion[[i]] = table(G_true, Estimated)
    confusion[[i]][1,1] = confusion[[i]][1,1]- p #non devo considerare la diagonale
    sbagliati = confusion[[i]][1,1] + confusion[[i]][2,2]
    sensitivity[i] = confusion[[i]][2,2] / (confusion[[i]][2,1] + confusion[[i]][2,2])
    if(confusion[[i]][1,1] == 0) {
      uno_meno_specificity[i] = 0
    } else {
      uno_meno_specificity[i] = confusion[[i]][1,1] / (confusion[[i]][1,1] + confusion[[i]][1,2])
    }
    giusti[i] =  sbagliati / (p*p - p)
  }

  best_soglia = tol[which.max(giusti)]
  best_cut    = which.max(giusti)
  best_confusion = confusion[[best_cut]]
  best_graph   = matrix(0,p,p)
  best_graph[PL >= best_soglia]   = 1

  if(isTRUE(ROC)) {
    #Plot ROC
    x11()
    plot(x = (1 - uno_meno_specificity), y = sensitivity, type = 'b',  col = 'red' , pch = 16, lwd = 2,
         main = 'ROC Curve', xlab = '1 - spec', ylab = 'sens')
    text((1 - uno_meno_specificity),  sensitivity, tol, col = 'black', cex = 0.6, pos = 4)
  }
  result = list()
  result[[1]] = confusion
  result[[2]] = best_soglia
  result[[3]] = best_confusion
  result[[4]] = best_graph
  result[[5]] = compute_AUC(uno_meno_specificity, sensitivity)
  names(result) = c('all_confusion', 'best_threshold', 'best_confusion', 'best_truncated_graph', 'AUC')
  return(result)

}

###########################################################################

###########################################################################
# hamming_distance_analysis -----------------------------------------------

#' Bayesian Hamming Distance Analysis
#'
#' Given the true graph and a plinks matrix, this utiliry computes the hamming distances between
#' the real graph and the truncated ones.
#' @param PL The posterior probability of inclusion of each link according to \code{\link{FGM_sampler}}.
#' It should be an upper diagonal matrix.
#' @param G_true The true graph.
#' @param tol The sequence of thresholds, used to select a graph truncating PL at that value.
#' @return A list of 3 elements:
#' \itemize{
#'   \item{\strong{best_threshold}, the best value of \code{tol} according to this analysis.}
#'   \item{\strong{best_truncated_graph}, the proposed posterior graph according to the analysis.}
#'   \item{\strong{best_distance}, the best value of the hamming distance, the one corresponding to \code{best_threshold}.}
#' }

#' @export
hamming_distance_analysis = function(PL, G_true, tol = seq(0.1, 1, by = 0.05)) {
  #suppressMessages(library(e1071))
  p = dim(PL)[1]
  G_true_vett = G_true[upper.tri(G_true, diag = F)]
  if(any(tol > max(PL)))
    tol <- tol[-which(tol > max(PL))]
  if(is.null(tol))
    stop("No feasible tolerances")
  distanze = rep(0,length(tol))
  for (i in 1:length(tol)) {
    tolerance <- tol[i]
    Estimated = matrix(0,p,p)
    Estimated[PL >= tolerance] = 1
    Estimated = Estimated[upper.tri(Estimated, diag = F)]
    distanze[i] =  e1071::hamming.distance(Estimated, G_true_vett)
  }
  best_soglia_hamm = tol[which.min(distanze)]
  best_cut_ham     = distanze[which.min(distanze)]
  best_graph_ham   = matrix(0,p,p)
  best_graph_ham[PL >= best_soglia_hamm]   = 1
  best_graph_ham   = best_graph_ham + t(best_graph_ham)
  best_soglia_hamm

  result = list()
  result[[1]] = best_soglia_hamm
  result[[2]] = best_graph_ham
  result[[3]] = best_cut_ham
  names(result) = c('best_treshold', 'best_truncated_graph', 'best_distance')
  return(result)
}

###########################################################################

###########################################################################
# FDR_analysis ------------------------------------------------------------

#' Bayesian FDR Analysis
#'
#' Given the plinks matrix, this utility computes the False Discovery Rate Index,
#' forcing the false discovery rates to be less than \code{0.05}-
#' @param PL The posterior probability of inclusion of each link according to \code{\link{FGM_sampler}}.
#' It should be an upper diagonal matrix.
#' @param tol The sequence of thresholds, used to select a graph truncating PL at that value.
#' @param min_rate The upper bound for the false discoveries.
#' @return a list of 3 elements:
#' \itemize{
#'   \item{\strong{best_threshold}, the best value of \code{tol} according to this analysis.}
#'   \item{\strong{best_truncated_graph}, the proposed posterior graph according to the analysis.}
#'   \item{\strong{best_distance}, the best FDR computed.}
#'}

#' @export
FDR_analysis = function(PL, tol = seq(0.01, max(PL), by = 0.01), min_rate = 0.05) {
  PL_vet = PL[upper.tri(PL, diag = F)]
  if(any(tol > max(PL)))
    tol <- tol[-which(tol > max(PL))]
  if(is.null(tol))
    stop("No feasible tolerances")

  FDR = rep(0,length(tol))

  for (i in 1:length(tol)) {
    tolerance <- tol[i]
    sopra_soglia = PL_vet[PL_vet >= tolerance]
    FDR[i] = sum( 1 - sopra_soglia )/length(sopra_soglia)
  }

  if(FDR[1] < min_rate){
    best_soglia_fdr = tol[1]
   }else
      for(i in 2:length(FDR)){
        if(FDR[i] < min_rate)
          break()
      }

  best_soglia_fdr = tol[i]
  FDR[i]
  best_soglia_fdr
  best_graph_fdr = matrix(0,dim(PL)[1],dim(PL)[2])
  best_graph_fdr[PL >= best_soglia_fdr]   = 1
  best_graph_fdr <- best_graph_fdr + t(best_graph_fdr)

  result = list()
  result[[1]] = best_soglia_fdr
  result[[2]] = best_graph_fdr
  names(result) = c('best_treshold', 'best_truncated_graph')
  return(result)

}

###########################################################################

###########################################################################
# threshold_analysis ------------------------------------------------------
#' Multiple a posteriori model selection
#'
#' This function performs in a single call the sensitivity, the hamming distance,and the FDR analysis.
#' See \code{\link{misclass_analysis}}, \code{\link{hamming_distance_analysis}} and \code{\link{FDR_analysis}} for information.
#' @usage threshold_analysis(PL, G_true, tol = seq(0.1, 1, by = 0.05),  min_rate = 0.05, ROC = FALSE)
#' @param PL The posterior probability of inclusion of each link according to \code{\link{FGM_sampler}}.
#' It should be an upper diagonal matrix.
#' @param G_true The true graph.
#' @param tol The sequence of thresholds, used to select a graph truncating PL at that value.
#' @param min_rate The upper bound for the false discoveries.
#' @param ROC Boolean parameter. If \code{TRUE}, the plot the ROC curve is showed.
#' @return A list with the outputs of sensitivity, hamming distance and FDR analysis.

#' @export
threshold_analysis <- function(PL, G_true, tol = seq(0.1, 1, by = 0.05),  min_rate = 0.05, ROC = FALSE) {
  result <- list()
  result[[1]] <- misclass_analysis(PL, G_true, tol, ROC)
  result[[2]] <- hamming_distance_analysis(PL, G_true, tol)
  result[[3]] <- FDR_analysis(PL, tol, min_rate)
  names(result) <- c('sensitivity', 'hamming', 'FDR')
  return(result)
}

###########################################################################
