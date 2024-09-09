### Graph of graphs: Concentration of the untransformed posterior

# Code for the plot on concentration of the untransformed posterior in the paper
# "Graph of Graphs: From Nodes to Supernodes in Graphical Models" by Maria De
# Iorio, Willem van den Boom, Alexandros Beskos, Ajay Jasra and Andrea
# Cremaschi.


source("graph_of_graphs_MCMC.R")
p <- 6L
n <- 60L
prec_mat <- diag(p)
node1 <- 1:(p %/% 2L)
node2 <- (p %/% 2L + 1L):p

for (i in node1[-1]) {
  prec_mat[i - 1L, i] <- 0.4
  prec_mat[i, i - 1L] <- prec_mat[i - 1L, i]
}

for (i in node2[-1]) {
  prec_mat[i - 1L, i] <- prec_mat[1, 2]
  prec_mat[i, i - 1L] <- prec_mat[i - 1L, i]
}

prec_mat[1, p %/% 2L + 1L] <- 0.0
prec_mat[p %/% 2L + 1L, 1] <- prec_mat[1, p %/% 2L + 1L]
set.seed(1L)
Y <- mvtnorm::rmvnorm(n = n, sigma = solve(prec_mat))
log_cohesion = dnbinom(x = 0:(p - 1L), size = 2, prob = 0.5, log = TRUE)


compute_exact_posterior <- function (Y) {
  # Compute the posterior by exhaustive enumeration.
  df_0 <- 3
  Y_cor <- abs(cor(Y))
  log_w_matrix <- compute_log_w_matrix(X = Y, coars_pow = 1)
  p <- NCOL(Y)
  n <- NROW(Y)
  
  memoise_mar_lik_supernode <- memoise::memoise(f = function(ind) {
    compute_log_tree_activation_f(
      log_w_matrix = log_w_matrix[ind, ind], p = length(ind)
    )
  }, cache = cachem::cache_mem(max_size = Inf))
  
  memoise_posterior_no_G <- memoise::memoise(
    f = function (centers, memoise_mar_lik_supernode) {
      compute_log_posterior_no_G(
        Y, Y_cor, centers, memoise_mar_lik_supernode, log_cohesion,
        nested = FALSE
      )
    },
    cache = cachem::cache_mem(max_size = Inf),
    omit_args = "memoise_mar_lik_supernode"
  )
  
  memoise_Laplace <- memoise::memoise(f = function (adj, U) log_p_GGM_Laplace(
    adj = adj, df_0 = 3, U = U, n = n
  ), cache = cachem::cache_mem(max_size = Inf))
  
  # All possible center configurations
  centers_mat <- as.matrix(expand.grid(replicate(
    n = p, expr = c(FALSE, TRUE), simplify = FALSE
  )))[-1, ]
  
  # All possible supergraphs
  G_list <- vector(mode = "list", length = p)
  
  for (n_centers in 1:p) {
    n_superedges <- (n_centers * (n_centers - 1L)) %/% 2L
    n_supergraphs <- as.integer(2^n_superedges)
    
    G_list[[n_centers]] <- array(
      data = 0L, dim = c(n_supergraphs, n_centers, n_centers)
    )
    
    if (n_superedges == 0L) next
    
    G_mat <- as.matrix(expand.grid(replicate(
      n = n_superedges, expr = 0:1, simplify = FALSE
    )))
    
    upper_ind <- upper.tri(G_list[[n_centers]][1, , ])
    
    for (G_ind in 1:n_supergraphs) {
      G_list[[n_centers]][G_ind, , ][upper_ind] <- G_mat[G_ind, ]
      
      G_list[[n_centers]][G_ind, , ] <- G_list[[n_centers]][G_ind, , ] +
        t(G_list[[n_centers]][G_ind, , ])
    }
  }
  
  n_C <- dim(centers_mat)[1]
  log_prob_centers <- numeric(n_C)
  log_pmf <- matrix(data = NA_real_, nrow = n_C, ncol = dim(G_list[[p]])[1])
  
  for (centers_ind in 1:n_C) {
    centers <- centers_mat[centers_ind, ]
    n_centers <- sum(centers)
    n_supergraphs <- dim(G_list[[n_centers]])[1]
    log_prob_G <- numeric(n_supergraphs)
    
    for (G_ind in 1:n_supergraphs) {
      G <- G_list[[n_centers]][G_ind, , ]
      
      tmp <- compute_log_posterior(
        centers, G, memoise_mar_lik_supernode, memoise_posterior_no_G,
        memoise_Laplace, coars_pow_lik = 1, nested = FALSE, pi_se = 0.5
      )
      
      log_prob_G[G_ind] <- tmp$log_ret
      assignment <- centers2tessellation(Y_cor, centers)
    }
    
    log_sum_prob_G <- matrixStats::logSumExp(log_prob_G)
    log_prob_centers[centers_ind] <- log_sum_prob_G
    log_pmf[centers_ind, 1:n_supergraphs] <- log_prob_G
  }
  
  # Normalize the posterior
  log_norm_const <- matrixStats::logSumExp(log_prob_centers)
  prob_centers <- exp(log_prob_centers - log_norm_const)
  log_pmf <- log_pmf - log_norm_const
  
  # Compute posterior coclustering probabilities and posterior on the number of
  # centers
  post_sim <- matrix(data = 0, nrow = p, ncol = p)
  
  for (centers_ind in 1:n_C) {
    centers <- centers_mat[centers_ind, ]
    tmp <- centers2tessellation(Y_cor, centers)
    
    for (i in 1:p) for (j in 1:p) if (tmp[i] == tmp[j]) {
      post_sim[i, j] <- post_sim[i, j] + prob_centers[centers_ind]
    }
    
    n_centers <- sum(centers)
  }
  
  return (list(
    post_sim = post_sim, log_pmf = log_pmf
  ))
}


n_vec <- 10L * 1:6
n_n <- length(n_vec)
res_arr <- array(data = NA_real_, dim = c(n_n, p, p))
compute_edge_prob <- TRUE

for (i in 1:n_n) {
  print(i)
  res_arr[i, , ] <- compute_exact_posterior(Y = Y[1:n_vec[i], ])$post_sim
}


plot_sim_mat <- function (mat = NULL) {
  p <- NCOL(mat)
  diag(mat) <- 1
  mat <- mat[, p:1]
  colnames(mat) <- p:1
  rownames(mat) <- 1:p
  n_levels <- 1000L
  
  lattice::levelplot(
    x = mat,
    col.regions = gray(n_levels:0 / n_levels),
    xlab = "Node",
    ylab = "Node",
    at = seq(from = 0, to = 1, length.out = n_levels + 1L),
    scales = list(draw = FALSE),
    colorkey = list(
      space = "bottom", title = "Posterior coclustering probability"
    ),
    panel = function (...) {
      lattice::panel.levelplot(...)
      
      lattice::panel.abline(h = p - 3L + 0.5, col = "red", lty = 2L)
      lattice::panel.abline(v = 3L + 0.5, col = "red", lty = 2L)
    }
  )
}


library(latticeExtra)
pdf("conc_sim_mat.pdf", width = 7, height = 2.7)

print(c(
  "n = 10" = plot_sim_mat(res_arr[1, , ]),
  "n = 20" = plot_sim_mat(res_arr[2, , ]),
  "n = 30" = plot_sim_mat(res_arr[3, , ]),
  "n = 40" = plot_sim_mat(res_arr[4, , ]),
  "n = 50" = plot_sim_mat(res_arr[5, , ]),
  "n = 60" = plot_sim_mat(res_arr[6, , ]),
  layout = c(n_n, 1L), merge.legends = FALSE
))

dev.off()
