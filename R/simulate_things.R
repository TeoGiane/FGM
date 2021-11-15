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

## Little library that collects utilities useful for simulations ##

###########################################################################
# FGM_sim -----------------------------------------------------------------

#' Functional Graphical Model simulation
#'
#' This utility simulates data from the Functional Graphical Model under analysis. Since the seed of
#' the random generator can be fixed, this allows multiple replications of the same simulated data
#' in order to perform multiple check on the performed analysis.
#'
#' @param data.seed Natural number, the seed for the random generator for the data.
#' @param graph.seed Natural number, the seed for the random generator for the Graph-Precision Matrix pair.
#' @param X.axis An r-dimensional vector that contains the abscissa points in which the data are evaluated.
#' @param dims A named vector collecting the dimensions of the input data, i.e. n, number of curves, r, number of gridpoints and p, the basis dimension.
#' @param basis A string that declares which type of basis needs to be used in order to embed data in their
#' functional framework. Right now, only "Spline" is available, ad uses a b-spline basis of third order.
#' @param sparsity \eqn{[0, 1]}-ranged parameter that indicates how sparse the precision matrix of the
#' functional graphical model is required. Default value is \eqn{0.3}.
#' @param type Select the desider Peterson graph/precision matrix couple. See \code{\link{KG_sim}}.
#' @return A list containing all the simulated data for the Functional Graphical Model in the form of an FGM object and a list with the parameters of simulation, that can be useful for debugging purposes.

#' @export
FGM_sim <- function(data.seed, graph.seed, X.axis, dims, basis = "Spline", method, sparsity = 0.3, type = 4) {

  if(is.null(data.seed) || is.null(graph.seed) || is.null(X.axis) || is.null(dims) || is.null(method)){
    stop("One or more inputs are missing, with no default.")
  }

  if(is.na(dims['n']) || is.na(dims['r']) || is.na(dims['p'])){
    stop("The dims vector is not well formed.")
  }

  # Unpacking dims vector
  n <- dims['n']
  r <- dims['r']
  p <- dims['p']

  # Generate basemat
  if (basis == "Spline") {
    m = 3; degree = m-1
    basis = fda::create.bspline.basis(rangeval = range(X.axis), nbasis=p, norder=m)
    basemat = fda::eval.basis(X.axis, basis)
  } else {
    stop("Unrecognized type of spline required.")
  }

  tau_eps <- 100
  set.seed(data.seed)
  mu <- as.numeric(FGM::rmvnorm(n = 1, mean = rep(0, p), sigma = 0.001*diag(p)))

  list_G_K <- KG_sim(method, p, sparsity, graph.seed, type)
  G <- list_G_K$G
  K <- list_G_K$K

  beta <- matrix(0, nrow = n, ncol = p)
  data <- matrix(0, nrow = n, ncol = r)
  K_inv <- solve(K)
  for (i in 1:n) {
    set.seed(data.seed + i)
    beta[i, ] <- FGM::rmvnorm(n = 1, mean = mu, sigma = K_inv)
    set.seed(data.seed + 10*i)
    data[i, ] <- FGM::rmvnorm(n = 1, mean = basemat %*% beta[i, ], sigma = diag(rep(1/tau_eps, r)))
  }

  # Creaing the param list
  params <- list()
  params[[1]] <- beta
  params[[2]] <- mu
  params[[3]] <- tau_eps
  params[[4]] <- G
  params[[5]] <- K
  names(params) <- c("beta", "mu", "tau_eps", "G", "K")

  # Creating the simulated_data list
  simulated_data <- list()
  simulated_data[[1]] <- FGM(X.axis, data, basis = "Spline", dim.basis = p)
  simulated_data[[2]] <- params
  names(simulated_data) <- c("FGMobj", "params")

  # Return function output
  return(simulated_data)

}

###########################################################################

###########################################################################
# KG_sim ------------------------------------------------------------

#' Graph and Precision matrix joint simulation
#'
#' This function provides a graph simulation at a fixed seed for a random generator, according
#' to different methods. In particular, three methods are here implemented:
#' \itemize{
#'   \item{\code{Bernoulli} associates at each link generation a Bernoulli trial.}
#'   \item{\code{Cluster} induces a graph with a blocking structure, according to the sparsity parameter.}
#'   \item{\code{Peterson} this method exploits literature precisions matrices commonly used for simulations.
#'   For further references, read, for instance, \url{https://www.researchgate.net/publication/261020849}.}
#' }
#' @param method The String identifying which method will be used to simulate the graph.
#' @param p The B-spline basis dimension. Is required to be multiple of 4, if \code{method="Peterson"}.
#' @param sparsity The [0, 1]-ranged parameter that indicates how sparse the precision matrix of the
#' functional graphical model is required. Used in case \code{method='Bernoulli'}, otherwise discarded.
#' @param S0 Natural number, the seed for the random generator for the data.
#' @param type A natural number between 1 and 4, which selects which Peterson matrix will be used
#' for the sampling. See \code{\link{Peterson_p16}} and \code{\link{Peterson_p40}}.
#' @return A list containing the sampled graph and the induced precision matrix.

#' @export
KG_sim <- function(method, p, sparsity, S0, type) {

  if(method == 'Bernoulli') {
    set.seed(S0)
    G <- BDgraph::bdgraph.sim(p = p, graph = "random", n = 1, prob = sparsity)$G
    G <- str2adj(p, adj2str(G))

    set.seed(S0)
    K <- BDgraph::rgwish(adj = G, b = 4)

    list_G_K <- list()
    list_G_K[[1]] <- G

    set.seed(S0)
    K <- BDgraph::rgwish(adj = G, b = 4)

    list_G_K[[2]] <- K
    names(list_G_K) <- c('G', 'K')

    return(list_G_K)

  }

  if(method == 'Cluster') {
    set.seed(S0)
    G <- bdgraph.sim(p = p, graph = "cluster", n = 1, prob = sparsity)$G
    G <- str2adj(p, adj2str(G))

    set.seed(S0)
    K <- BDgraph::rgwish(adj = G, b = 4)

    list_G_K <- list()
    list_G_K[[1]] <- G

    set.seed(S0)
    K <- rgwish(adj = G, b = 4)

    list_G_K[[2]] <- K
    names(list_G_K) <- c('G', 'K')

    return(list_G_K)
  }

  if(method == 'Peterson') {
    if( p%%4 != 0)
      stop('p should be multiple of 4')
    string <- paste0('Peterson_p',p)
    data(list = string); list_G_K <- get(string);
    if(type == 1)
      return(list_G_K$`1`)
    if(type == 2)
      return(list_G_K$`2`)
    if(type == 3)
      return(list_G_K$`3`)
    if(type == 4)
      return(list_G_K$`4`)
  } else
    stop("Method does not match")
}

###########################################################################
