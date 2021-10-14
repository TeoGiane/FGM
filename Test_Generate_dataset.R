# Set seeds for the data generation
Nsims <- 1
graph.seed <- 96
set.seed(graph.seed)
data.seed <- extraDistr::rdunif(Nsims, graph.seed^2, graph.seed^3)


# Bernoulli 03 ------------------------------------------------------------

# Generate datasets and save
p <- 7
n <- 200
data('purees')
X = purees$wavelengths
range_x <- range(X)
r = length(X)
dims = c("n"=n, "r"=r, "p"=p)
for (i in 1:Nsims) {
  dati_simulati <- FGM_sim(data.seed=data.seed[i], graph.seed=graph.seed,
                           X.axis = X, dims = dims,
                           method = "Bernoulli", sparsity = 0.3)
  filename <- sprintf("Test_p7_%02d.Rdat", i)
  save(dati_simulati, file = filename)
}
