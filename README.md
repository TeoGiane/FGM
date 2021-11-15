# FGM - Functional Graphical Models
**Authors:** Laura Codazzi, Alessandro Colombi, Matteo Gianella.

The `R` package **FGM** provides statistical tools for Bayesian structural learning in the framework of functional data analysis. The package implements the model presented in
[Codazzi et al. (2021)](https://arxiv.org/abs/2103.11666). To learn the structre of dependencies, this package exploits the efficien Birth and Death approach presented in  
[Mohammadi and Wit (2015)](https://projecteuclid.org/journals/bayesian-analysis/volume-10/issue-1/Bayesian-Structure-Learning-in-Sparse-Gaussian-Graphical-Models/10.1214/14-BA889.full)
and implemented in the `R` package `BDgraph`, available [here](https://github.com/cran/BDgraph). Such a package is not loaded in **FGM** as the authors needed to slightly modify
the **BDgraph** code to adapt it to their needs. As a consequence, parts of the `BDgraph` code have been included in the package. For those, all rights are reserved to `BDgraph`.

## Installation
**FGM** requires the loading of some `R` packages. Run in the `R` console the following line to make sure that they are available.
```R
install.packages(c("stringr", "coda", "fields", "e1071", "fda"))
```
The package can be installed using:
```R
devtools::install_github("TeoGiane/FGM")
library("FGM")
```

## Example
`FGM_sampler` is the function where the sampling strategy is implemented. Here is an example of its usage to fit the Purees dataset available [here](https://data.mendeley.com/datasets/frrv2yd9rg/1).
We considered only the subgroup made of strawberry purees, discarding the others.
```R
# load data ---------------------------------------------------------------

data("StrawberryPurees")
data("StrawberryWavelengths")
fgmobj <- create_structure()

# set parameters ----------------------------------------------------------

#dim
n = dim(StrawberryPurees)[1]
r = dim(StrawberryPurees)[2]
p = 40

#prior
a_Beta = 1
b_Beta = (2*p - 2)/3 - 1 #Beta prior, Beta(a_Beta, b_Beta)
gprior <- c(a_Beta, b_Beta) # Graph prior

#hyperparameters
hyper <- list()
hyper[[1]] <- 2
hyper[[2]] <- 0.02
hyper[[3]] <- 100
hyper[[4]] <- 3
hyper[[5]] <- gprior
names(hyper) <- c('a_tau_eps', 'b_tau_eps', 'sigma_mu', 'd0', 'gprior')

# initial values ----------------------------------------------------------

init <- list()
init[[1]] = fgmobj$basemat
init[[2]] = diag(rep(1,fgmobj$p))
init[[3]] = fgmobj$beta
init[[4]] = fgmobj$mu
init[[5]] = diag(rep(1,fgmobj$p))
init[[6]] = 1
names(init) = c('basemat', 'G0', 'Beta0', 'mu0', 'K0', 'tau_eps0')

# sampler parameters ------------------------------------------------------

niter <- 100000
nburn <-  20000
thin  <-     20
thinG <-      2
iter_to_store <- (niter - nburn)/thin
iter_to_storeG <- (niter - nburn)/thinG

# run ---------------------------------------------------------------------

fit <- FGM_sampler(data = data, niter=niter, nburn=nburn, thin=thin, thinG=thinG, init=init, hyper=hyper)

# read results ------------------------------------------------------------

beta_list = list()
for(i in 1:dim(fgmobj$basemat)[1]){
  beta_list[[i]] = fit[[i]]
}

beta_mean = beta_pointwise_estimate(beta_list, p)
mu_mean = colMeans(fit$mu)
tau_eps_mean = mean(fit$tau_eps)
PL = plinks(fit$bdobject)
G_fdr = FDR_analysis(PL, tol = seq(0.1,1,by = 0.01))$best_truncated_graph
K_hat = fit$bdobject$K_hat

```
