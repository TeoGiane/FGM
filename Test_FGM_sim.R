# Load data ---------------------------------------------------------------
filename <- "Test_p7_01.Rdat"
load(filename)

# Extract problem dimension -----------------------------------------------
n <- as.numeric(dati_simulati$FGMobj$dims['n'])
r <- as.numeric(dati_simulati$FGMobj$dims['r'])
p <- as.numeric(dati_simulati$FGMobj$dims['p'])

# Setting prior parameters ------------------------------------------------
gprior <- c(0.5,0.5)
hyper <- list()
hyper[[1]] <- 2 * 10
hyper[[2]] <- 2 * 0.001
hyper[[3]] <- 100
hyper[[4]] <- 3
hyper[[5]] <- matrix(0, p, p)
hyper[[6]] <- gprior
names(hyper) <- c('a_tau_eps', 'b_tau_eps', 'sigma_mu', 'd0', 'D0', 'gprior')
init = dati_simulati$FGMobj$init

# Setting sampler parameter -----------------------------------------------
niter <- 6000
nburn <- 1000
thin  <- 20
thinG <- 1
iter_to_store <- (niter - nburn)/thin
iter_to_storeG <- (niter - nburn)/thinG

# Sampler execution -------------------------------------------------------
fit <- FGM_sampler(data = dati_simulati$FGMobj$data,
                   niter=niter, nburn=nburn, thin=thin, thinG=thinG,
                   init=init, hyper=hyper)



# see results -------------------------------------------------------------

PL = plinks(fit$bdobject)
library(fields)
x11();image.plot(PL)
PL
dati_simulati$params$G
