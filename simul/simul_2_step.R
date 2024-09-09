### Graph of graphs: Comparison with two-step approaches

# Code for the simulation studies with two-step approaches in the paper "Graph
# of Graphs: From Nodes to Supernodes in Graphical Models" by Maria De Iorio,
# Willem van den Boom, Alexandros Beskos, Ajay Jasra and Andrea Cremaschi.



# Install required packages
pkgs <- c("igraph", "CVglasso", "latticeExtra", "mvtnorm")
pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(pkgs) > 0L) {
  print("Installing required packages...")
  install.packages(pkgs = pkgs, repos = "https://cloud.r-project.org")
}



get_both <- function (X) {
  # Two-step approach
  p <- ncol(X)
  set.seed(1L)
  print("Fitting glasso...")
  
  adj_matrix <- CVglasso::CVglasso(
    X = X, lam.min.ratio = 1e-5,
    # Five is the default value for the number K of folds for cross-validation.
    cores = 5L
  )$Omega != 0
  
  # Since `glasso` might produce an assymmetric `Omega`, we only use the upper
  # triangle of `Omega`/`adj_matrix` as suggested in Rolfs and Rajaratnam
  # (2013, doi:10.1016/j.csda.2012.07.013).
  adj_matrix[lower.tri(adj_matrix)] <- t(adj_matrix)[lower.tri(adj_matrix)]
  
  print("Clustering nodes...")
  
  graph <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj_matrix, mode = "undirected", diag = FALSE
  )
  
  mat_betweenness <- matrix(0, nrow = p, ncol = p)
  
  for (ind_vec in igraph::communities(
    igraph::cluster_edge_betweenness(graph)
  )) mat_betweenness[ind_vec, ind_vec] <- 1
  
  mat_eigen <- matrix(0, nrow = p, ncol = p)
  
  for (ind_vec in igraph::communities(igraph::cluster_leading_eigen(graph))) {
    mat_eigen[ind_vec, ind_vec] <- 1
  }
  
  return (list(
    adj_matrix = adj_matrix, mat_betweenness = mat_betweenness,
    mat_eigen = mat_eigen
  ))
}


get_lat_fac_data <- function (n, supernode_sizes, superedge, vareps) {
  p <- sum(supernode_sizes)
  K <- length(supernode_sizes)
  
  # Latent factors from multivariate Gaussian
  Psi <- diag(K)
  if (superedge) {
    # Introduce dependence between the latent factors.
    Psi[1,3] <- Psi[3,1] <- 0.9
  }
  
  lat_fac <- mvtnorm::rmvnorm(n = n, mean = rep(0,K), sigma = solve(Psi))  
  
  X <- matrix(NA, nrow = n, ncol = p)
  tmp <- c(1L, 1L + cumsum(supernode_sizes))
  
  for (k in 1:K) for (j in tmp[k]:(tmp[k + 1] - 1L)) {
    X[, j] <- lat_fac[, k] + rnorm(n = n, sd = sqrt(vareps))
  }
  
  return (x = scale(X, center = TRUE))
}



supernode_sizes <- c(5L, 10L, 20L)
burnin <- 15000
n_iter <- 5000

n_simulations <- 4
res_list <- vector("list", length = n_simulations)




print("Latent factors")

print("Without superedge")

set.seed(123)
print("n = 1000, vareps = 0.01")
n <- 1000

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = FALSE, vareps = 0.01
)

res_list[[1]] <- get_both(X)

set.seed(123)
print("n = 1000, vareps = 0.05")
n <- 1000

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = FALSE, vareps = 0.05
)

res_list[[2]] <- get_both(X)

print("With superedge")

set.seed(123)
print("n = 1000, vareps = 0.01")
n <- 1000

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = TRUE, vareps = 0.01
)

res_list[[3]] <- get_both(X)

set.seed(123)
print("n = 1000, vareps = 0.05")
n <- 1000

X <- get_lat_fac_data(
  n = n, supernode_sizes = supernode_sizes, superedge = TRUE, vareps = 0.05
)

res_list[[4]] <- get_both(X)


plot_mat <- function (
    mat = NULL
) {
  p <- NCOL(mat)
  mat <- mat[, p:1]
  
  lattice::levelplot(
    x = mat,
    col.regions = c("white", "black"),
    xlab = "Node",
    ylab = "Node",
    scales = list(draw = FALSE),
    colorkey = FALSE,
    panel = function (...) {
      lattice::panel.levelplot(...)
      
      for (ind in cumsum(c(5L, 10L, if (p == 75L) 20L else NULL))) {
        lattice::panel.abline(h = p - ind + 0.5, col = "red", lty = 2L)
        lattice::panel.abline(v = ind + 0.5, col = "red", lty = 2L)
      }
    }
  )
}


library(latticeExtra)

pdf(
  file = paste0("simul_lat_fac_2_step.pdf"),
  width = 10.5, height = 14.5
)

print(c(
  "With superedge, var = 0.05: Graph" = plot_mat(res_list[[4]]$adj_matrix),
  "Edge betweenness" = plot_mat(res_list[[4]]$mat_betweenness),
  "Leading eigenvector" = plot_mat(res_list[[4]]$mat_eigen),
  "With superedge, var = 0.05: Graph" = plot_mat(res_list[[3]]$adj_matrix),
  "Betweennes" = plot_mat(res_list[[3]]$mat_betweenness),
  "Leading eigenvector" = plot_mat(res_list[[3]]$mat_eigen),
  "Without superedge, var = 0.01: Graph" = plot_mat(res_list[[2]]$adj_matrix),
  "Edge betweenness" = plot_mat(res_list[[2]]$mat_betweenness),
  "Leading eigenvector" = plot_mat(res_list[[2]]$mat_eigen),
  "Without superedge, var = 0.01: Graph" = plot_mat(res_list[[1]]$adj_matrix),
  "Edge betweenness" = plot_mat(res_list[[1]]$mat_betweenness),
  "Leading eigenvector" = plot_mat(res_list[[1]]$mat_eigen),
  layout = c(3L, 4L)
))

dev.off()
