# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C)  Laura Codazzi, Alessandro Colombi, Matteo Gianella                           |
#                                                                                                 |
#     FGM is free software: you can redistribute it and/or modify it under                        |
#     the terms of the GNU General Public License as published by the Free                        |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
#                                                                                                 |
#     Maintainer: Laura Codazzi (laura.codazzi@tuhh.de),                                          |
#                 Alessandro Colombi (a.colombi10@campus.unimib.it),                              |
#                 Matteo Gianella (matteo.gianella@polimi.it)                                     |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

## Little library to manipulate the beta list ouput of extract_all_beta ##

###########################################################################
# beta_pointwise_estimate -------------------------------------------------
#' Compute pointwise estimate of beta
#'
#' This function provides the pointwise estimate of any beta w.r.t. the usual
#' quadratic loss function.
#'
#' @param beta The list of beta coming from select_all_beta function or a
#' dataframe corresponding to a specific curve.
#' @param p The number of basis.
#'
#' @return An nxp matrix containing the pointwise estimate of each beta or
#' a p vector in case a single beta dataframe is passed as input.

#' @export
beta_pointwise_estimate = function(beta, p){
  if(is.data.frame(beta)){
    result = colMeans(as.matrix(beta))
    names(result) <- NULL
    return(result)
  }
  else if(is.list(beta)) {
    n = length(beta)
    result <- matrix(0, nrow = n, ncol = p)
    for (i in 1:n) {
      beta_media <- colMeans(beta[[i]])
      result[i,] <- beta_media
    }
    return(result)
  }
  else {
    stop('Wrong input inserted.')
  }
}

###########################################################################

###########################################################################
# SS_beta -----------------------------------------------------------------
#' Compute the sum of squares of beta
#'
#' This utility simply compute the sum of sqares of the betas related to the
#' i-th curve. Useful tool in out custom GS.
#'
#' @param beta The list of beta coming from select_all_beta function.
#' @param mu The p-dimensional vector of the mean.
#'
#' @return A p-dimensional vector a for which a[j] = \ifelse{html}{\out{&Sigma;<sub>i</sub>(&beta;<sub>ij</sub> - &mu;<sub>j</sub>)<sup>2</sup>}}{\eqn{\sum_{i}\left(\beta_{ij}-\mu_{j}\right)^2}}.

#' @export
SS_beta = function(beta, mu = NULL){
  n = length(beta)
  p = dim(beta[[1]])[2]
  if(is.null(mu)){
    mu = rep(0,p)
  }
  beta_matrix = as.matrix(beta[[1]]); colnames(beta_matrix) <- NULL
  for(i in 2:n){
    tmp_mat <- as.matrix(beta[[i]]); colnames(tmp_mat) <- NULL
    beta_matrix = rbind(beta_matrix, tmp_mat)
  }
  result = rep(0, p)
  for(j in 1:p){
    result[j] = sum( (beta_matrix[ ,j] - mu[j])^2)
  }
  return(result)
}

###########################################################################

###########################################################################
# mean_beta_per_component -------------------------------------------------
#' Compute the component-wise mean for beta
#'
#' This function compute the mean value of beta per component.
#'
#' @param beta The list of beta coming from select_all_beta function.
#'
#' @return A p-dimensional vector containing the mean value of the
#' j-th component of beta.

#' @export
mean_beta_per_component = function(beta){
  n = length(beta)
  beta_matrix = as.matrix(beta[[1]]); colnames(beta_matrix) <- NULL
  for(i in 2:n){
    tmp_mat <- as.matrix(beta[[i]]); colnames(tmp_mat) <- NULL
    beta_matrix = rbind(beta_matrix, tmp_mat)
  }
  beta_medie = colMeans(beta_matrix)
  return(beta_medie)
}

###########################################################################
