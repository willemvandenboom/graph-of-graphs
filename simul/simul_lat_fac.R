### Graph of graphs: Simulation based on latent factors

# Code for the simulation studies based on latent factors in the paper
# "Graph of Graphs: From Nodes to Supernodes in Graphical Models" by Maria De
# Iorio, Willem van den Boom, Alexandros Beskos, Ajay Jasra and Andrea
# Cremaschi.


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
    pi_se_prop = pi_se_prop
  )
  
  coars_pow <- 1 / n
  
  print("Running nested MCMC")
  
  res_nested <- run_MCMC(
    X = X, nested = TRUE, n_iter = n_iter, burnin = burnin, n_thin = n_thin,
    log_cohesion = log_cohesion, coars_pow = coars_pow, pi_se = pi_se,
    pi_se_prop = pi_se_prop
  )
  
  return (list(
    post_sim_joint = res_joint$post_sim,
    superedge_mat_joint = res_joint$superedge_mat,
    post_sim_nested = res_nested$post_sim,
    superedge_mat_nested = res_nested$superedge_mat
  ))
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

n_simulations <- 4
res_list <- vector("list", length = n_simulations)


print("Without superedge")

set.seed(123)
print("vareps = 0.01")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = FALSE, vareps = 0.01
)

res_list[[1]] <- get_both(X)

set.seed(123)
print("vareps = 0.05")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = FALSE, vareps = 0.05
)

res_list[[2]] <- get_both(X)

print("With superedge")

set.seed(123)
print("vareps = 0.01")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = TRUE, vareps = 0.01
)

res_list[[3]] <- get_both(X)

set.seed(123)
print("vareps = 0.05")

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = TRUE, vareps = 0.05
)

res_list[[4]] <- get_both(X)


plot_sim_mat <- function (mat) {
  p <- NCOL(mat)
  mat <- mat[, p:1]
  n_levels <- 1000L
  
  lattice::levelplot(
    x = mat,
    col.regions = gray(n_levels:0 / n_levels),
    xlab = "Node",
    ylab = "Node",
    at = seq(from = 0, to = 1, length.out = n_levels + 1L),
    scales = list(draw = FALSE),
    par.strip.text = list(cex = 0.85),
    colorkey = if (label == "Posterior coclustering probability") {
      list(space = "bottom", title = label)
    } else {
      FALSE
    },
    panel = function (...) {
      lattice::panel.levelplot(...)
      
      for (ind in cumsum(c(5L, 10L, if (p == 75L) 20L else NULL))) {
        lattice::panel.abline(h = p - ind + 0.5, col = "red", lty = 2L)
        lattice::panel.abline(v = ind + 0.5, col = "red", lty = 2L)
      }
    }
  )
}


pdf_width <- 7
pdf_height <- 14.5
library(latticeExtra)


label <- "Posterior coclustering probability"
pdf(file = "simul_lat_fac_sim_mat.pdf", width = pdf_width, height = pdf_height)
print(c(
  "Superedge, n = 1000, var = 0.05: Coars. lik." = plot_sim_mat(res_list[[4]]$post_sim_joint),
  "Superedge, n = 1000, var = 0.05: Nested MCMC" = plot_sim_mat(res_list[[4]]$post_sim_nested),
  "Superedge, n = 1000, var = 0.01: Coars. lik." = plot_sim_mat(res_list[[3]]$post_sim_joint),
  "Superedge, n = 1000, var = 0.01: Nested MCMC" = plot_sim_mat(res_list[[3]]$post_sim_nested),
  "No superedge, n = 1000, var = 0.05: Coars. lik." = plot_sim_mat(res_list[[2]]$post_sim_joint),
  "No superedge, n = 1000, var = 0.05: Nested MCMC" = plot_sim_mat(res_list[[2]]$post_sim_nested),
  "No superedge, n = 1000, var = 0.01: Coars. lik." = plot_sim_mat(res_list[[1]]$post_sim_joint),
  "No superedge, n = 1000, var = 0.01: Nested MCMC" = plot_sim_mat(res_list[[1]]$post_sim_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()


label <- "Supergraph"

pdf(
  file = "simul_lat_fac_superedge_medianG.pdf",
  width = pdf_width, height = pdf_height
)

print(c(
  "Superedge, n = 1000, var = 0.05: Coars. lik." = plot_sim_mat(round(res_list[[4]]$superedge_mat_joint)),
  "Superedge, n = 1000, var = 0.05: Nested MCMC" = plot_sim_mat(round(res_list[[4]]$superedge_mat_nested)),
  "Superedge, n = 1000, var = 0.01: Coars. lik." = plot_sim_mat(round(res_list[[3]]$superedge_mat_joint)),
  "Superedge, n = 1000, var = 0.01: Nested MCMC" = plot_sim_mat(round(res_list[[3]]$superedge_mat_nested)),
  "No superedge, n = 1000, var = 0.05: Coars. lik." = plot_sim_mat(round(res_list[[2]]$superedge_mat_joint)),
  "No superedge, n = 1000, var = 0.05: Nested MCMC" = plot_sim_mat(round(res_list[[2]]$superedge_mat_nested)),
  "No superedge, n = 1000, var = 0.01: Coars. lik." = plot_sim_mat(round(res_list[[1]]$superedge_mat_joint)),
  "No superedge, n = 1000, var = 0.01: Nested MCMC" = plot_sim_mat(round(res_list[[1]]$superedge_mat_nested)),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()
