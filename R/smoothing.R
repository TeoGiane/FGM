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

## This library contains the utility to perform a bayesian smoothing ##

#' Bayesian smoothing for the Functional (Graphical) Model
#'
#' This function computes the smoothed curves and their credibility bands, if required.
#' @param beta A list of length \eqn{n}, with all the sampled values for each beta coefficient. One can either
#' pass the post-processed output for the Functional Model coming from \code{\link{select_all_beta}} or
#' the first \eqn{n} elements of the output list for the Functional Graphical Model, coming from \code{\link{FGM_sampler}}.
#' @param basemat The design matrix given by the evaluation of the basis over the grid of points.
#' @param alpha A [0-1]-ranged parameter to be set in order to get the \ifelse{html}{\out{(1-&alpha;)\%}}{\eqn{(1-\alpha)\%}} credibility bands.
#' @param band Boolean, \code{FALSE} by default. If \code{TRUE}, the function also computes the \ifelse{html}{\out{(1-&alpha;)\%}}{\eqn{(1-\alpha)\%}}
#' credibility bands.
#' @return In case band is \code{FALSE}, it returns a matrix with all the values of the smoothed curves.
#' If band is \code{TRUE}, it returns a list of three matrices collecting the lower, the mean and the upper bands values of the smoothing.

#' @export
smooth_FGM = function(beta, basemat, alpha = 0.1, band = FALSE){
  if(!is.list(beta)){
    stop('Beta should be a list with all the sampled values.')
  }
  niter = dim(beta[[1]])[1]
  p = dim(beta[[1]])[2]
  n = length(beta)
  t_mean = beta_pointwise_estimate(beta, p)
  r = dim(basemat)[1]

  if(band == TRUE){
    result = list()
    lower_band = matrix(0, nrow = n, ncol = r)
    upper_band = matrix(0, nrow = n, ncol = r)
    mean = matrix(0, nrow = n, ncol = r)
    for(i in 1:n){
      tmp <- sapply(beta[[i]], function(x) quantile(x, c(alpha/2,1-alpha/2)))
      y_band <- sapply(as.data.frame(t(tmp)), function(x) basemat %*% x)
      y_mean = basemat%*%t_mean[i,]
      lower_band[i, ] = y_band[,1]
      upper_band[i, ] = y_band[,2]
      mean[i, ] = y_mean
    }
    result[[1]] = lower_band
    result[[2]] = upper_band
    result[[3]] = mean
    names(result) = c('lower_band', 'upper_band', 'mean')
    return(result)
  }else{
    result = matrix(0, nrow = n, ncol = r)
    for(i in 1:n){
      y_mean = basemat%*%t_mean[i,]
      result[i,] = y_mean
    }
    return(result)
  }
}
