### Graph of graphs: Comparison with a similarity function representing a
### complete graph

# Code for the simulation studies with a similarity function representing a
# complete graph in the paper "Graph of Graphs: From Nodes to Supernodes in
# Graphical Models" by Maria De Iorio, Willem van den Boom, Alexandros Beskos,
# Ajay Jasra and Andrea Cremaschi.


source("graph_of_graphs_MCMC.R")


get_both <- function (X) {
  # Run MCMC with coarsening of the likelihood and nested MCMC.
  p <- ncol(X)
  n <- nrow(X)
  log_cohesion = dnbinom(x = 0:(p-1), size = 2, prob = 1 / p, log = TRUE)
  coars_pow <- 10 / n
  pi_se <- 0.1
  pi_se_prop <- 0.1
  
  print("Running MCMC with coarsening of the likelihood")
  
  res_joint <- run_MCMC(
    X = X, nested = FALSE, n_iter = n_iter, burnin = burnin, n_thin = n_thin,
    log_cohesion = log_cohesion, coars_pow = coars_pow, pi_se = pi_se,
    pi_se_prop = pi_se_prop, use_t = TRUE
  )
  
  print("Posterior on the number of supernodes:")
  print(table(rowSums(res_joint$centers_MCMC)) / n_iter)
  
  coars_pow <- 1 / n
  
  print("Running nested MCMC")
  
  res_nested <- run_MCMC(
    X = X, nested = TRUE, n_iter = n_iter, burnin = burnin, n_thin = n_thin,
    log_cohesion = log_cohesion, coars_pow = coars_pow, pi_se = pi_se,
    pi_se_prop = pi_se_prop, use_t = TRUE
  )
  
  print("Posterior on the number of supernodes:")
  print(table(rowSums(res_nested$centers_MCMC)) / n_iter)
}


get_lat_fac_data <- function (n, supernode_sizes, superedge, vareps) {
  p <- sum(supernode_sizes)
  K <- length(supernode_sizes)
  
  # Latent factors from multivariate Gaussian
  Psi <- diag(K)
  
  if (superedge) {
    # Introduce dependence between the latent factors.
    Psi[1, 3] <- Psi[3, 1] <- 0.9
  }
  
  lat_fac <- mvtnorm::rmvnorm(n = n, mean = rep(0, K), sigma = solve(Psi))  
  X <- matrix(NA, nrow = n, ncol = p)
  tmp <- c(1L, 1L + cumsum(supernode_sizes))
  
  for (k in 1:K) for (j in tmp[k]:(tmp[k + 1] - 1L)) {
    X[, j] <- lat_fac[, k] + rnorm(n = n, sd = sqrt(vareps))
  }
  
  return (x = scale(X, center = TRUE))
}



supernode_sizes <- c(5L, 10L, 20L)
n <- 1000L
burnin <- 15000
n_iter <- 5000
n_thin <- 5


print("Without superedge")

set.seed(123)
print("vareps = 0.01")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = FALSE, vareps = 0.01
)

get_both(X)

set.seed(123)
print("vareps = 0.05")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = FALSE, vareps = 0.05
)

get_both(X)

print("With superedge")

set.seed(123)
print("vareps = 0.01")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = TRUE, vareps = 0.01
)

get_both(X)

set.seed(123)
print("vareps = 0.05")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = TRUE, vareps = 0.05
)

get_both(X)
