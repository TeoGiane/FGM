## This library collect the Gibbs Sampler for our Functional Graphical Model ##

#' Sampling from the posterior distribution of the Functional Graphical Model
#'
#' This utility implements an hybrid Gibbs Sampler strategy in order to sample
#' the posterior distribution of the parameters \ifelse{html}{\out{<em>&beta;<sub>i</sub>, &mu;, K, G, &tau;<sub>&epsilon;</sub></em>}}{\eqn{\{\beta_{i}\}_{i}, \mu, K, G, \tau_{\epsilon}}}
#' of the Functional Graphical Model for the spectrometric data analysis performed.
#'
#' @param data The data that needs to be fitted through this model.
#' @param niter Number of iterations for the sampler.
#' @param nburn Number of burn-in iterations.
#' @param thin Thining parameter.
#' @param thinG Thining parameter for the joint precision-graph distribution.
#' @param init List that collects the starting point of the Gibbs Sampler Chains.
#' \itemize{
#'   \item{\strong{basemat}, the B-spline rxp basis matrix.}
#'   \item{\strong{G0}, starting graph (See \code{\link{bdgraph}} function for more information about the \code{g.start} input).}
#'   \item{\strong{Beta0}, list of p-dimensional vector, starting point for the \eqn{\beta} chains.}
#'   \item{\strong{mu0}, p-dimensional vector, starting point of the \eqn{\mu} chain.}
#'   \item{\strong{K0}, pxp Precision matrix, starting point for the \eqn{K} chain.}
#'   \item{\strong{tau_eps0}, real, starting point for the \ifelse{html}{\out{<em>&tau;<sub>&epsilon;</sub></em>}}{\eqn{\tau_\epsilon}} chain.}
#' }
#' @param hyper List that collects all the needed hyperparameters for the model. If nothing is passed, hyperparameters are set at their default values.
#' \itemize{
#'   \item{\strong{a_tau_eps}, first parameter of gamma prior for \ifelse{html}{\out{<em>&tau;<sub>&epsilon;</sub></em>}}{\eqn{\tau_\epsilon}}. \cr Default value: \eqn{2*10}.}
#'   \item{\strong{b_tau_eps}, second parameter of gamma prior for \ifelse{html}{\out{<em>&tau;<sub>&epsilon;</sub></em>}}{\eqn{\tau_\epsilon}}. \cr Default value: \eqn{2*0.001}.}
#'   \item{\strong{sigma_mu}, mean for the Gaussian prior for \eqn{\mu}. \cr Default value: \eqn{100}.}
#'   \item{\strong{d0}, degrees of freedom of the G-Wishart prior for the precision matrix. \cr Default value: \eqn{3}.}
#'   \item{\strong{gprior}, prior probability of each link of the graph (See \code{\link{bdgraph}} function for more information about the \code{g.prior} input). \cr Default value: \eqn{0.5}.}
#' }
#' @return An \eqn{n+4}-item list with the following components:
#' \itemize{
#'   \item{\strong{beta_i} with \eqn{i=1:n}, the first n components are \code{(iter_to_store x p)} matrices storing the draws from the posterior distribution for all the \ifelse{html}{\out{<em>&beta;<sub>i</sub></em>}}{\eqn{\{\beta_{i}\}_{i}}}.}
#'   \item{\strong{mu}, an \code{(iter_to_store x p)} matrix collecting the \ifelse{html}{\out{<em>&mu;</em>}}{\eqn{\mu}} draws from the posterior distribution.}
#'   \item{\strong{tau_eps}, an \code{(iter_to_store)} vector collecting the \ifelse{html}{\out{<em>&tau;<sub>&epsilon;</sub></em>}}{\eqn{\tau_{\epsilon}}} draws from the posterior distribution.}
#'   \item{\strong{bdobject}, a \code{bdgraph} object (see \code{\link{bdgraph}} for details).}
#'   \item{\strong{K}, an \code{(iter_to_storeG)} vector collecting the \eqn{K} draws from the posterior distribution.}
#' }

#' @export
FGM_sampler <- function(data, niter, nburn, thin, thinG, init, hyper=NULL) {
  #library('BDgraph')
  if( niter < 0 )
  stop('niter must be positive')
  # Check data
  if(is.null(data) || is.null(niter) || is.null(nburn) ||
     is.null(thin) || is.null(thinG) || is.null(init)) {
    stop("One or more parameters are missing, with no default.")
  }

  # Unpacking the init list
  if(is.list(init)){
    basemat <- init$basemat
    G0 <- init$G0
    Beta0 <- init$Beta0
    mu0 <- init$mu0
    K0 <- init$K0
    tau_eps0 <- init$tau_eps0
  } else {
    stop("No list is provided as init parameter.")
  }

  # Setting dimensions
  p <- dim(basemat)[2]
  r <- dim(basemat)[1]
  n <- dim(data)[1]

  # Check if custom hyperparameters are provided
  if(is.null(hyper)) {
    a_tau_eps <- 2*10
    b_tau_eps <- 2*0.001
    sigma_mu <- 100
    d0 <- 3
    gprior <- 0.5
  } else if (is.list(hyper)) {
    # Unpacking the hyper list
    a_tau_eps <- hyper$a_tau_eps
    b_tau_eps <- hyper$b_tau_eps
    sigma_mu <- hyper$sigma_mu
    d0 <- hyper$d0
    gprior <- hyper$gprior
  } else {
    stop("No list is provided as custom hyper parameter.")
  }

  # Building the result structure, built in this way:
  # result[[i]] for i=1:n contains the (iter_to_store x p)-matrix with the Beta_i draws
  # result[[n+1]] contains the (iter_to_store x p) matrix with the mu draws
  # result[[n+2]] contains a vector with the tau_eps draws
  # result[[n+3]] contains a bdgraph object, storing the information on K and G
  # result[[n+4]] contains a vector with the K draws, in a string form
  iter_to_store = (niter - nburn)/thin
  iter_to_storeG = (niter - nburn)/thinG
  M = matrix(0, nrow = iter_to_store, ncol = p)
  X = data.frame(row.names = '1', codice = ' ', peso = 0)
  result = list()
  result2 <- rep("", iter_to_storeG)
  #result2 = list()
  for(i in 1:n){
    result[[i]] = M                       # beta_i, for i=1:n
  }
  result[[n+1]] = M                       # mu
  result[[n+2]] = rep(0,iter_to_store)    # tau_eps

  # Adding bdgraph construction helpers in the list
  result[[n+3]] = X                       # sampled_graphs and weights
  result[[n+4]] = rep(0, iter_to_storeG)  # all visited graphs
  result[[n+5]] = matrix(0, p, p)         # K_hat
  #for(i in 1 : iter_to_storeG){
  #  result2[[i]] = matrix(0, p, p)
  #}
  all_weights = rep(0, iter_to_storeG)

  # Fixed quantities
  tbase_base = t(basemat) %*% basemat   # Used to update betas
  base_y = list()                       # base_y[[i]] = t(basemat)*y_i, used to update betas as well
  Si = rep(0,n)                         # Si[i] = t(y_i) %*% y_i, used to update tau_eps
  for(i in 1:n){
    #base_y[[i]] = t(basemat)%*%as.numeric(data[i,])
    #Si[i] = sum(as.numeric(data[i,])^2)
    base_y[[i]] = t(basemat) %*% as.numeric(data[i,])
    Si[i] = sum(data[i,]^2)
  }

  # Fixed updated hyperparameters
  a_eps <- (n*r + a_tau_eps)/2

  # Chain initialization
  beta <- Beta0
  mu <- mu0
  K <- K0
  tau_eps <- tau_eps0

  cat('\nFGM - Hybrid Gibbs Sampler started:\n')
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  #cat('Start')
  #cat('\n')
  for(iter in 1:niter){

    # beta_i parameters update and draw, for i=1:n
    B_n <- solve(tau_eps*tbase_base + K)
    Kmu = K%*%mu
    for (i in 1:n) {
      bn_i <- B_n %*% (tau_eps*base_y[[i]] + Kmu)
      beta[[i]] <- FGM::rmvnorm(1, mean = bn_i, sigma = B_n)
    }


    # mu parameters update and draw
    Sni = rep(0,p)
    for(i in 1:n){
      Sni = Sni + beta[[i]]
    }
    A = solve(diag(rep(1/sigma_mu, p)) + n*K)
    #inv_A = solve(A)
    #a = inv_A %*% (K%*%t(Sni))
    #mu = as.numeric(FGM::rmvnorm(1, mean = a, sigma = inv_A))
    a = A %*% (K%*%t(Sni))
    mu = as.numeric(FGM::rmvnorm(1, mean = a, sigma = A))

    # (K, G) parameters update and joint draw
    data_BD = matrix(0,nrow = p, ncol = p)
    for(i in 1:n){
      data_BD = data_BD  + t(beta[[i]] - mu)%*%(beta[[i]] - mu)
    }
    if(iter != 1) {
      fit = FGM::bdgraph(data = data_BD, n = n, method = 'ggm', algorithm = 'bdmcmc', iter = 1, burnin = 0,
                    g.prior = gprior, df.prior = d0, g.start = fit$last_graph, save = T)
    } else {
      fit = FGM::bdgraph(data = data_BD, n = n, method = 'ggm', algorithm = 'bdmcmc', iter = 1, burnin = 0,
                    g.prior = gprior, df.prior = d0, g.start = G0, save = T)
    }
    K = fit$last_K

    # tau_eps parameters update and draw
    temp <- 0
    for(i in 1:n) {
      #y_i <- as.numeric(data[i, ])
      y_i <- data[i, ]
      #temp <- temp +  t( y_i - basemat %*% t(beta[[i]]) ) %*% (y_i - basemat %*% t(beta[[i]]))  #cosi devo fare piu conti
      temp <- temp +  Si[i] +  beta[[i]]%*%tbase_base%*%t(beta[[i]]) -2 * beta[[i]]%*%base_y[[i]]
    }
    c_eps <- 0.5 * (b_tau_eps + temp)
    tau_eps <- rgamma(1, a_eps, c_eps)
    #print('ciao')

    # Storing result
    # beta, mu, tau_eps storage
    if(iter > nburn & (iter - nburn) %% thin==0) {
      stored = (iter - nburn) / thin
      # beta
      for(i in 1:n){
        result[[i]][stored, ] = beta[[i]]
      }
      #mu
      result[[n + 1]][stored, ] = mu
      # tau_eps
      result[[n + 2]][stored] = tau_eps
    }

    # (K,G) storage
    if(iter > nburn & (iter - nburn) %% thinG==0) {
      stored = (iter - nburn) / thinG
      # (K,G)
      if(stored == 1) {
        result[[n+3]]$codice = adj2str(fit$last_graph)
        result[[n+3]]$peso   = fit$all_weights
        result[[n+4]][1] = 1
      } else {
        nuovo.grafo = adj2str(fit$last_graph)
        if(!any( result[[n+3]]$codice ==  nuovo.grafo )){ # If TRUE, the visited graph is new
          nuova.riga = data.frame(row.names = as.character(dim(result[[n+3]])[1] + 1 ),
                                  codice = nuovo.grafo, peso = fit$all_weights)
          result[[n+3]] = rbind(result[[n+3]], nuova.riga)
          result[[n+4]][stored] = dim(result[[n+3]])[1]
        } else {
          index = which(result[[n+3]]$codice == nuovo.grafo )
          result[[n+3]]$peso[index] = result[[n+3]]$peso[index] + fit$all_weights
          result[[n+4]][stored] = index
        }
      }
      all_weights[stored] = fit$all_weights
      result[[n + 5]] = result[[n+5]] + fit$K_hat*as.numeric(fit$graph_weights)
      #result2[[stored]] = toString(K[upper.tri(K,diag = T)])
      result2[stored] = toString(K[upper.tri(K,diag = T)])
    }
    setTxtProgressBar(pb, iter)
  }
  result[[n+5]] = result[[n+5]]/sum(result[[n+3]]$peso)

  # bdgraph object creation through helpers
  bdobj = FGM::bdgraph(diag(rep(4,5)), n = 3, iter = 1, save = T)
  bdobj$sample_graphs = result[[n+3]]$codice
  bdobj$graph_weights = result[[n+3]]$peso
  bdobj$K_hat = result[[n+5]]
  bdobj$all_graphs = result[[n+4]]
  bdobj$all_weights = all_weights
  bdobj$last_graph = str2adj(p,result[[n+3]]$codice[result[[n+4]][iter_to_store]])
  #bdobj$last_K = result2[[iter_to_storeG]]
  bdobj$last_K = str2K(p,result2[iter_to_storeG])

  # Cleaning the output
  result[[n+3]] = NULL
  result[[n+4]] = NULL
  result[[n+5]] = NULL

  # Adding last elements to the list
  result[[n+3]] = bdobj
  result[[n+4]] = result2

  # Setting names to the output list
  names(result) <- c(paste0("beta_",1:n), "mu", "tau_eps", "bdobject", "K")

  close(pb)
  return(result)
}
