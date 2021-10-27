## Utility that re-arrange data for the Gibbs Sampler ##

###########################################################################
# FGM ---------------------------------------------------------------------
#' FGM object constructor
#'
#' This utility recieves in input functional data defined in the format (X,Y), where X are the absciasse
#' at which the functions are evaluated and Y are the values of the functions in X. For a clarifying example,
#' see the format in which the \code{\link{purees}} dataset is stored in this package. It returns an FGM
#' object that collects useful information and helps the usage of these data in \code{\link{FGM_sampler}}.
#'
#' @param X.axis An r-dimensional vector that contains the abscissa points in which the data are evaluated.
#' @param Y.axis This can be both an r-dimensional vector f(X) or an n x p data frame which stores, in each
#' row, the values \ifelse{html}{\out{Y<sub>i</sub>=f(X<sub>i</sub>), i=1:n}}{\eqn{Y_i = f(X_i),~i=1:n}}.
#' @param basis A string that declares which type of basis needs to be used in order to embed data in their
#' functional framework. Right now, only "Spline" is available, ad uses a b-spline basis of third order.
#' @param dim.basis Integer, the dimension of the spline basis.
#'
#' @return A list of the elements, i.e.
#' \itemize{
#'   \item{\strong{data}, the dataframe provided in Y.axis.}
#'   \item{\strong{dims}, a named vector collecting the dimensions of the input data, i.e. n, number of curves, r, number of gridpoints and p, the basis dimension.}
#'   \item{\strong{init}, the \code{init} list that can be used as input in \code{\link{FGM_sampler}}.}
#' }
#'

#' @export
FGM <- function(X.axis, Y.axis, basis = "Spline", dim.basis) {

  # Set dimensions
  if(is.data.frame(Y.axis) || is.matrix(Y.axis)){
    n <- dim(Y.axis)[1]
    r <- dim(Y.axis)[2]
  } else {
    n <- 1
    r <- length(Y.axis)
  }
  p <- dim.basis
  dims <- c(n,r,p); names(dims) <- c("n", "r", "p")

  # Create Spline Basis
  if (basis == "Spline") {
    m = 3; degree = m-1
    basis = fda::create.bspline.basis(rangeval = range(X.axis), nbasis=p, norder=m)
    basemat = fda::eval.basis(X.axis, basis)
  } else {
    stop("Unrecognized type of spline required.")
  }

  # Generate auxiliary sturctures
  aux_vec <- rep(0,p)
  aux_matrix <- diag(rep(1,p))
  aux_list <- list()
  for (i in 1:n) {
    aux_list[[i]] <- aux_vec
  }

  # Create init list
  init <- list()
  init[[1]] <- basemat
  init[[2]] <- BDgraph::bdgraph.sim(p = p, graph = "random")$G
  init[[3]] <- aux_list
  init[[4]] <- aux_vec
  init[[5]] <- aux_matrix
  init[[6]] <- 1
  names(init) <- c("basemat", "G0", "Beta0", "mu0", "K0", "tau_eps0")

  # Create FGM output
  result <- list()
  result[[1]] <- Y.axis
  result[[2]] <- dims
  result[[3]] <- init
  names(result) <- c("data", "dims", "init")

  # Return
  return(result)
}

###########################################################################

###########################################################################
# create_structure --------------------------------------------------------
#' Data preprocessing for the FGM Hybrid Gibbs Sampler
#'
#' This function pre-process the data stored in the \code{purees} data frame
#' (see \code{\link{purees}} for more information) and creates a data structure
#' which will be useful for the sampling algorithm (i.e. \code{\link{FGM_sampler}}).
#'
#' @param n Number of curves which will be analysed. Default value is \code{351}.
#' @param p Dimension of the B-spline basis. Default value is \code{40}
#' @param wave.indexes Wavelengths range to analyse, expressed in terms of index.
#' Default value is the whole range of \code{235} wavelengths.
#'
#' @return A list with the re-arranged parameters provided in input.

#' @export
create_structure <- function(n = 351, p = 40, wave.indexes = c(1:235)){
  #library("stringr")
  #library("fda")

  result <- list()
  result[[1]] <- n
  result[[2]] <- p

  data('purees');data_tot = purees
  X <- data_tot$wavelengths
  strawberry <- data_tot$data[which(data_tot$data$Group=="Strawberry"),]
  Y <- as.matrix(strawberry[,-1])
  data <- Y[1:n, ]
  result[[3]] <- data


  range_x <- range(X)
  r = length(X)
  result[[4]] <- r
  m = 3; degree = m-1         # Spline Order m and Degree m-1
  basis = fda::create.bspline.basis(rangeval = range_x, nbasis=p, norder=m)
  basemat = fda::eval.basis(X, basis)

  result[[5]] <- basemat

  temp <- list()
  vec <- rep(0, p)
  for (i in 1:n) {
    temp[[i]] <- vec
  }
  result[[6]] <- temp
  result[[7]] <- vec
  result[[8]] <- 0

  names(result) <- c('n', 'p', 'data', 'r', 'basemat', 'beta', 'mu', 'tau')
  return(result)
}

###########################################################################
