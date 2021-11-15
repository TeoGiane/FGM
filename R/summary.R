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

## Useful functions for quick visualization of posterior results ##

###########################################################################
# summary_FM ----------------------------------------------------------
#' Summarizing posterior distribution of the FM model
#'
#' This utility provides useful plots for an overall summary of the posterior
#' distribution of all parameters which come from the JAGS output of the Functional Model
#'
#' @param output The JAGS output.
#' @param n.curves The number of curves for which we've computed the
#' posterior distribution.
#' @param alpha Real number used to compute the (1-\eqn{\alpha}) credible intervals.
#' Default value is 0.05, hence we compute the 95\% CI.
#'
#' @return A total of five plots which represent, in order, the mean
#' value of each beta in the form of a matrix, the correlation
#' plot of the beta coefficients, the (1-\eqn{\alpha})\% credible intervals
#' for the mu vector, the tau vector and the tau_eps parameter.

#' @export
summary_FM <- function(output, n.curves, alpha=0.05) {

  # Preparing result list
  result <- list()

  # Plot for mean value of each beta coefficient
  cat("Summarizing information about beta coefficients...\n")
  beta_list <- select_all_beta(output, n.curves)
  beta_means <- unlist(lapply(beta_list, function(x) colMeans(x)))
  beta_means <- matrix(beta_means, nrow = n.curves, ncol = length(beta_means)/n.curves, byrow = TRUE)
  p <- length(beta_means)/n.curves
  fields::image.plot(x=1:p, y=1:n.curves, t(beta_means),
             xlab="Bands", ylab='Curves', graphics.reset = TRUE)
  title(main="Mean Values for beta coefficients")
  result[[1]] <- beta_means

  # Pause between plots
  readline(prompt="Press [enter] to continue ")

  # Plot of credible intervals for mu
  cat("Summarizing information about mu coefficients...\n")
  mu_sample <- select_posterior_FM(output, "mu")
  mu_quant <- apply(mu_sample, 2, function(x) quantile(x, c(alpha/2,0.5,1-alpha/2)))
  plot(x=1:p, rep(NA,p), ylim = range(mu_quant), xlab='Index', ylab='')
  segments(x0=1:p, y0=mu_quant[1,], y1=mu_quant[3,], lwd=2, col='darkgrey')
  points(x=1:p,y=mu_quant[2,], pch=21, ylim = range(mu_quant), xlab='Index', ylab = '', col='black',bg='black')
  title(main=paste0((1-alpha)*100,'% Credible Intervals for mu'))
  result[[2]] <- mu_quant

  # Pause between plots
  readline(prompt="Press [enter] to continue ")

  # Plot of credible intervals for tau
  cat("Summarizing information about tau coefficients...\n")
  tau_sample <- select_posterior_FM(output, "tau")
  tau_quant <- apply(tau_sample, 2, function(x) quantile(x, c(alpha/2,0.5,1-alpha/2)))
  plot(x=1:p, rep(NA,p), ylim = range(tau_quant), xlab='Index', ylab='')
  segments(x0=1:p, y0=tau_quant[1,], y1=tau_quant[3,], lwd=2, col='darkgrey')
  points(x=1:p,y=tau_quant[2,], pch=21, ylim = range(mu_quant), xlab='Index', ylab = '', col='black',bg='black')
  title(main=paste0((1-alpha)*100,'% Credible Intervals for the precision vector tau'))
  result[[3]] <- tau_quant

  # Pause between plots
  readline(prompt="Press [enter] to continue ")

  # Plot of credible intervals for tau_eps
  cat("Summarizing information about tau_eps coefficient...\n")
  tau_sample <- select_posterior_FM(output, "tau_eps")
  tau_quant <- apply(tau_sample, 2, function(x) quantile(x, c(alpha/2,0.5,1-alpha/2)))
  plot(density(as.numeric(unlist(tau_sample))), lwd=2, col='blue', main='')
  abline(v=tau_quant[1], lwd=2, lty=2, col='darkgrey')
  abline(v=tau_quant[2], lwd=2)
  abline(v=tau_quant[3], lwd=2, lty=2, col='darkgrey')
  title(main=paste0((1-alpha)*100,'% Credible Intervals for tau_eps'))
  result[[4]] <- tau_quant

  # Set names
  names(result) <- c("beta", "mu", "tau", "tau_eps")
  return(result)
}

###########################################################################

###########################################################################
# summary_FGM ----------------------------------------------------------
#' Summarizing posterior distribution of parameters
#'
#' This utility provides useful plots for an overall summary of the posterior
#' distribution of all parameters which come from the our own Gibbs Sampler
#' called \code{Gibbs_model11}
#'
#' @param output The \code{\link{FGM_sampler}} output.
#' @param n.curves The number of curves for which we've computed the
#' posterior distribution.
#' @param alpha Real number used to compute the 1-alpha credible intervals.
#' Default value is 0.05, hence we compute the 95\% CI.
#'
#' @return A total of five plots which represent, in order, the mean
#' value of each beta in the form of a matrix, the correlation
#' plot of the beta coefficients, the \eqn{(1-alpha)\%} credible intervals
#' for the mu vector, the tau vector and the tau_eps parameter.

#' @export
summary_FGM <- function(output, n.curves, alpha=0.05) {

  # Preparing result list
  result <- list()

  # Plot for mean value of each beta coefficient
  cat("Summarizing information about beta coefficients...\n")
  p = dim(output[[1]])[2]
  beta_means <- matrix(0, nrow = n.curves, ncol = p, byrow = TRUE)
  for(i in 1:n.curves) {
    beta_means[i, ] = beta_pointwise_estimate(as.data.frame(output[[i]]), p)
  }
  fields::image.plot(x=1:p, y=1:n.curves, t(beta_means),
             xlab="Bands", ylab='Curves', graphics.reset = TRUE)
  title(main="Mean Values for beta coefficients")
  result[[1]] <- beta_means

  # Pause between plots
  readline(prompt="Press [enter] to continue")

  # Plot of credible intervals for mu
  cat("Summarizing information about mu coefficients...\n")
  #output[[n.curves + 1]], qui c'è mu, evito di copiarlo due volte
  mu_quant <- apply(output[[n.curves + 1]], 2, function(x) quantile(x, c(alpha/2,0.5,1-alpha/2)))
  plot(x=1:p, rep(NA,p), ylim = range(mu_quant), xlab='Index', ylab='')
  segments(x0=1:p, y0=mu_quant[1,], y1=mu_quant[3,], lwd=2, col='darkgrey')
  points(x=1:p,y=mu_quant[2,], pch=21, ylim = range(mu_quant), col='black', bg='black')
  title(main=paste0((1-alpha)*100,'% Credible Intervals for mu'))
  result[[2]] <- mu_quant

  # Pause between plots
  readline(prompt="Press [enter] to continue")

  # Plot of credible intervals for tau_eps
  cat("Summarizing information about tau_eps coefficient...\n")
  tau_quant <- quantile(output[[n.curves + 2]], c(alpha/2,0.5,1-alpha/2))
  plot(density(as.numeric(unlist(output[[n.curves + 2]]))), lwd=2, col='blue', main='')
  abline(v=tau_quant[1], lwd=2, lty=2, col='darkgrey')
  abline(v=tau_quant[2], lwd=2)
  abline(v=tau_quant[3], lwd=2, lty=2, col='darkgrey')
  title(main=paste0((1-alpha)*100,'% Credible Intervals for tau_eps'))
  result[[3]] <- tau_quant

  # Pause between plots
  readline(prompt="Press [enter] to continue")

  # Graphical part
  cat("Summarizing information about the graphical part...\n")
  bdobj = fit$bdobject
  S = summary(bdobj)
  plotcoda(bdobj)
  PL = plinks(bdobj)
  p = dim(PL)[1]
  Median_graph = matrix(0,p,p)
  Median_graph[PL>0.5] = 1
  result[[4]] = PL
  result[[5]] = Median_graph
  result[[6]] = FDR_analysis(PL)

  # Setting names for the output
  names(result) = c("beta", "mu", "tau_eps", 'plinks', 'Median_graph', 'FDR_analysis')
  return(result)

}

###########################################################################






