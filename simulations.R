### Graph sphere: Simulation studies

# Code for the simulation studies in Appendix H of the paper "Graph Sphere: From
# Nodes to Supernodes in Graphical Models" by Willem van den Boom, Maria De
# Iorio, Alexandros Beskos and Ajay Jasra.



source("graph_sphere_MCMC.R")


get_both <- function (X) {
  # Run MCMC with coarsening of the likelihood and nested MCMC
  p <- ncol(X)
  log_cohesion = dnbinom(x = 0:(p - 1L), size = 2, prob = 1 / p, log = TRUE)
  
  set.seed(1L)
  print("Running MCMC with coarsening of the likelihood...")
  
  res_joint <- run_MCMC(
    X = X, nested = FALSE, n_iter = n_iter, log_cohesion = log_cohesion,
    prob_informed = 0.1
  )
  print("Running nested MCMC...")
  
  res_nested <- run_MCMC(
    X = X, nested = TRUE, n_iter = n_iter, log_cohesion = log_cohesion
  )
  
  return (list(
    post_sim_joint = res_joint$post_sim,
    superedge_mat_joint = res_joint$superedge_mat,
    post_sim_nested = res_nested$post_sim,
    superedge_mat_nested = res_nested$superedge_mat
  ))
}


get_lat_fac_data <- function (supernode_sizes, superedge) {
  p <- sum(supernode_sizes)
  K <- length(supernode_sizes)
  set.seed(1L)
  lat_fac <- matrix(data = rnorm(n * K), nrow = n, ncol = K)  # Latent factors
  
  if (superedge) {
    # Introduce dependence between the latent factors.
    lat_fac[, K] <- sqrt(0.05) * lat_fac[, K - 1L] + sqrt(0.95) * lat_fac[, K]
  }
  
  X <- matrix(data = NA_real_, nrow = n, ncol = p)
  tmp <- c(1L, 1L + cumsum(supernode_sizes))
  
  for (k in 1:K) for (j in tmp[k]:(tmp[k + 1] - 1L)) {
    X[, j] <- lat_fac[, k] + rnorm(n = n, sd = 0.01)
  }
  
  return (scale(x = X, center = FALSE))
}


get_cor_data <- function (supernode_sizes, cor_vec) {
  p <- sum(supernode_sizes)
  cor_mat <- diag(p)
  K <- length(supernode_sizes)
  center_ind <- c(1L, 1L + cumsum(supernode_sizes[-K]))
  
  for (k in 1:K) {
    ind <- center_ind[k] - 1L + 1:supernode_sizes[k]
    cor_mat[ind, ind] <- cor_vec[k]
  }
  
  diag(cor_mat) <- 1L
  i1 <- center_ind[K - 1L]
  i2 <- center_ind[K]
  cor_mat[i1, i2] <- -.3
  cor_mat[i2, i1] <- cor_mat[i1, i2]
  set.seed(1L)
  return (mvtnorm::rmvnorm(n = n, sigma = cor_mat))
}


get_prec_data <- function (superedge) {
  p <- sum(supernode_sizes)
  prec_mat <- diag(p)
  K <- length(supernode_sizes)
  center_ind <- c(1L, 1L + cumsum(supernode_sizes[-K]))
  
  # Supernodes consist of star trees
  for (k in 1:K) for (i in 1:(supernode_sizes[k] - 1L)) {
    i1 <- center_ind[k]
    i2 <- center_ind[k] + i
    
    # Absolute partial correlation is bounded by
    # 1 / sqrt(supernode_sizes[k] - 1L) in this setting of star trees.
    prec_mat[i1, i2] <- 0.9 / sqrt(supernode_sizes[k] - 1L)
    
    prec_mat[i2, i1] <- prec_mat[i1, i2]
  }
  
  if (superedge) {
    i1 <- center_ind[K - 1L]
    i2 <- center_ind[K]
    prec_mat[i1, i2] <- -.15
    prec_mat[i2, i1] <- prec_mat[i1, i2]
  }
  
  set.seed(1L)
  
  return (scale(
    x = mvtnorm::rmvnorm(n = n, sigma = solve(prec_mat)), center = FALSE
  ))
}


n <- 1e3L
supernode_sizes <- c(5L, 10L, 20L, 40L)
n_iter <- 1e4L
res_list <- list()


print("Latent factors")

print("Without superedge")
print("n = 1000...")
X <- get_lat_fac_data(supernode_sizes = supernode_sizes, superedge = FALSE)
res_list[[1]] <- get_both(X)
print("n = 100...")
res_list[[2]] <- get_both(X[1:100, ])
print("With superedge")
print("n = 1000...")
X <- get_lat_fac_data(supernode_sizes = supernode_sizes, superedge = TRUE)
res_list[[3]] <- get_both(X)
print("n = 100...")
res_list[[4]] <- get_both(X[1:100, ])


print("Correlation matrix")

print("p = 75")
print("High correlation")
print("n = 1000...")
X <- get_cor_data(supernode_sizes, c(0.9, 0.8, 0.5, 0.4))
res_list[[5]] <- get_both(X)
print("n = 100...")
res_list[[6]] <- get_both(X[1:100, ])
print("Low correlation")
res_list[[7]] <- get_both(get_cor_data(supernode_sizes, c(0.5, 0.4, 0.3, 0.2)))
print("p = 35. High correlation. n = 1000...")
res_list[[8]] <- get_both(get_cor_data(supernode_sizes[-4], c(0.9, 0.8, 0.5)))


print("Precision matrix")
supernode_sizes <- c(5L, 10L, 20L)

print("Without superedge")
print("n = 1000...")
X <- get_prec_data(superedge = FALSE)
res_list[[9]] <- get_both(X)
print("n = 100...")
res_list[[10]] <- get_both(X[1:100, ])
print("With superedge")
print("n = 1000...")
X <- get_prec_data(superedge = TRUE)
res_list[[11]] <- get_both(X)
print("n = 100...")
res_list[[12]] <- get_both(X[1:100, ])


plot_sim_mat <- function (
    mat = NULL
) {
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
    colorkey = list(space = "bottom", title = label),
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
pdf("simul_lat_fac_sim_mat.pdf", width = pdf_width, height = pdf_height)

print(c(
  "With superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[4]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[4]]$post_sim_nested),
  "With superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[3]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[3]]$post_sim_nested),
  "Without superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[2]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[2]]$post_sim_nested),
  "Without superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[1]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[1]]$post_sim_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()
pdf("simul_cor_sim_mat.pdf", width = pdf_width, height = pdf_height)

print(c(
  "High cor., n = 1000, p = 35: Coars. lik." = plot_sim_mat(res_list[[8]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[8]]$post_sim_nested),
  "Low cor., n = 1000, p = 75: Coars. lik." = plot_sim_mat(res_list[[7]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[7]]$post_sim_nested),
  "High cor., n = 100, p = 75: Coars. lik." = plot_sim_mat(res_list[[6]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[6]]$post_sim_nested),
  "High cor., n = 1000, p = 75: Coars. lik." = plot_sim_mat(res_list[[5]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[5]]$post_sim_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()
pdf("simul_prec_sim_mat.pdf", width = pdf_width, height = pdf_height)

print(c(
  "With superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[12]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[12]]$post_sim_nested),
  "With superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[11]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[11]]$post_sim_nested),
  "Without superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[10]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[10]]$post_sim_nested),
  "Without superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[9]]$post_sim_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[9]]$post_sim_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()


label <- "Posterior superedge probability"
pdf("simul_lat_fac_superedge.pdf", width = pdf_width, height = pdf_height)

print(c(
  "With superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[4]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[4]]$superedge_mat_nested),
  "With superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[3]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[3]]$superedge_mat_nested),
  "Without superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[2]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[2]]$superedge_mat_nested),
  "Without superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[1]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[1]]$superedge_mat_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()
pdf("simul_cor_superedge.pdf", width = pdf_width, height = pdf_height)

print(c(
  "High cor., n = 1000, p = 35: Coars. lik." = plot_sim_mat(res_list[[8]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[8]]$superedge_mat_nested),
  "Low cor., n = 1000, p = 75: Coars. lik." = plot_sim_mat(res_list[[7]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[7]]$superedge_mat_nested),
  "High cor., n = 100, p = 75: Coars. lik." = plot_sim_mat(res_list[[6]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[6]]$superedge_mat_nested),
  "High cor., n = 1000, p = 75: Coars. lik." = plot_sim_mat(res_list[[5]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[5]]$superedge_mat_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()
pdf("simul_prec_superedge.pdf", width = pdf_width, height = pdf_height)

print(c(
  "With superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[12]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[12]]$superedge_mat_nested),
  "With superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[11]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[11]]$superedge_mat_nested),
  "Without superedge, n = 100: Coars. lik." = plot_sim_mat(res_list[[10]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[10]]$superedge_mat_nested),
  "Without superedge, n = 1000: Coars. lik." = plot_sim_mat(res_list[[9]]$superedge_mat_joint),
  "Nested MCMC" = plot_sim_mat(res_list[[9]]$superedge_mat_nested),
  layout = c(2L, 4L), merge.legends = FALSE
))

dev.off()
