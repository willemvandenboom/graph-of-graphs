### Graph of graphs: Markov chain Monte Carlo

# Implementation of the Markov chain Monte Carlo described in the paper "Graph
# of Graphs: From Nodes to Supernodes in Graphical Models" by Maria De Iorio,
# Willem van den Boom, Alexandros Beskos, Ajay Jasra and Andrea Cremaschi.



# Install required packages
pkgs <- c(
  "abind", "BH", "CholWishart", "latticeExtra", "matrixStats", "memoise",
  "mvtnorm", "RcppBlaze", "withr"
)

pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(pkgs) > 0L) {
  print("Installing required packages...")
  install.packages(pkgs = pkgs, repos = "https://cloud.r-project.org")
}

df_0 <- 3  # Degrees of freedom of the G-Wishart prior for the supergraph



print("Compiling Rcpp code...")

# The following uses the Rcpp function `update_G_Rcpp` to provide an MCMC update
# for a graph in a Gaussian graphical model. This implements WWA
# (van den Boom et al., doi:10.1080/10618600.2022.2050250).
withr::with_makevars(
  new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
  code = Rcpp::sourceCpp(file = "src/ggm.cpp")
)


#' MCMC update for a Gaussian graphical model
#' 
#' @description The MCMC update is a marginal update on the graph in the
#'  Gaussian graphical model with a G-Wishart prior. The rate matrix of the
#'  G-Wishart prior is chosen as the identity matrix. The edges are assumed to
#'  be independent a-priori with equal prior inclusion probability.
#' 
#' @param adj The symmetric adjacency matrix of the graph with zeroes on the
#'   diagonal.
#' @param edge_prob The prior edge inclusion probability.
#' @param df_0 The degrees of freedom of the G-Wishart prior.
#' @param U The scatter matrix of the observations.
#' @param n The number of observations constituting `U`.
#' @return The symmetric adjacency matrix of the graph after the MCMC update.
update_G_WWA <- function (adj, edge_prob, df_0, U, n, n_edge) {
  p <- nrow(adj)
  
  return (tryCatch(
    expr = update_G_Rcpp(
      adj = adj,
      
      edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
      df = df_0 + n,
      df_0 = df_0,
      rate = diag(p) + U,
      n_edge = n_edge,  # Number of single edge updates attempted in this MCMC update
      seed = sample.int(n = .Machine$integer.max, size = 1L),
      loc_bal = FALSE
    ), error = function (e) {
      warning(paste("`update_G_Rcpp` failed:", e))
      return (adj)
    }
  ))
}


# Modified version of base::La.svd
fast_svd <- function (x, nu = min(n, p), nv = min(n, p), n = NULL, p = NULL) {
  if (is.null(n)) n <- nrow(x)
  if (is.null(p)) p <- ncol(x)
  
  if (nu || nv) {
    np <- min(n, p)
    if (nu <= np && nv <= np) {
      jobu <- "S"
      u <- matrix(0, n, np)
      vt <- matrix(0, np, p)
      nu0 <- nv0 <- np
    }
    else {
      jobu <- "A"
      u <- matrix(0, n, n)
      vt <- matrix(0, p, p)
      nu0 <- n
      nv0 <- p
    }
  }
  else {
    jobu <- "N"
    u <- matrix(0, 1L, 1L)
    vt <- matrix(0, 1L, 1L)
  }
  res <- .Internal(La_svd(jobu, x, double(min(n, p)), u, vt))
  res <- res[c("d", if (nu) "u", if (nv) "vt")]
  if (nu && nu < nu0) 
    res$u <- res$u[, seq_len(min(n, nu)), drop = FALSE]
  if (nv && nv < nv0) 
    res$vt <- res$vt[seq_len(min(p, nv)), , drop = FALSE]
  res
}


# Compute all pairwise weights with coarsening.
compute_log_w_matrix <- function (X, coars_pow) {
  n <- NROW(X)  # Number of observations
  p <- NCOL(X)  # Number of nodes
  delta <- 3  # Degrees of freedom of the Hyper Inverse Wishart prior
  delta_star <- 3 + n
  D <- diag(p)
  D_star <- D + crossprod(X)
  
  compute_log_f <- function (delta, Dij) {
    lgamma(0.5 * (delta + 1)) + 0.5 * delta * (
      log(Dij[1, 1]) + log(Dij[2, 2])
    ) - lgamma(0.5 * delta) - 0.5 * (delta + 1) * determinant(Dij)$modulus
  }
  
  log_w_matrix <- matrix(data = NA_real_, nrow = p, ncol = p)
  
  for (i in 1:p) {
    # Diagonal elements are included to compute the normalizing constant.
    log_w_matrix[i, i] <- lgamma(0.5 * delta_star) +
      0.5 * delta * log(D[i, i]) - 0.5 * n * log(pi) - lgamma(0.5 * delta) -
      0.5 * delta_star * log(D_star[i, i])
    
    if (i < p) for (j in (i + 1L):p) {
      log_w_matrix[i, j] <- compute_log_f(
        delta_star, D_star[c(i, j), c(i, j)]
      ) - compute_log_f(delta, D[c(i, j), c(i, j)])
      
      log_w_matrix[j, i] <- log_w_matrix[i, j]
    }
  }
  
  return (coars_pow * log_w_matrix)
}


# Compute the tree activation function.
compute_log_tree_activation_f <- function (
    log_w_matrix, p, eps = 1e-12
) {
  if (p == 1L) log_diag_Laplacian <- -Inf else {
    log_diag_Laplacian <- numeric(p)
    
    for (i in 1:p) {
      log_diag_Laplacian[i] <- matrixStats::logSumExp(log_w_matrix[i, -i])
    }
  }
  
  # Compute the normalizing constant.
  log_ret <- (2L - p) * log(p) + sum(diag(as.matrix(log_w_matrix)))
  if (p == 2L) log_ret <- log_ret + log_diag_Laplacian[1]
  
  if (p > 2L) {
    # Version of the Laplacian with scaling of rows and columns to ensure unit
    # diagonal to reduce over- and underflow errors.
    Laplacian <- diag(p)
    
    log_det <- NULL
    
    for (i in 1:p) {
      Laplacian[i, -i] <- -exp(log_w_matrix[i, -i] - 0.5 * (
        log_diag_Laplacian[i] + log_diag_Laplacian[-i]
      ))
      
      # Reduce condition number of Laplacian by clipping its smallest values.
      ind <- which(-Laplacian[i, -i] < eps)
      n_ind <- length(ind)
      if (n_ind == 0L) next
      Laplacian[i, i] <- 1 + sum(Laplacian[i, ind]) + n_ind * eps
      Laplacian[i, ind] <- -eps
    }
    
    for (u in sample(1:p)) {
      svd_res <- fast_svd(
        x = Laplacian[-u, -u], nu = 0L, nv = 0L, n = p - 1L, p = p - 1L
      )
      
      if (all(svd_res$d > 0)) {
        log_det <- sum(log(svd_res$d)) + sum(log_diag_Laplacian[-u])
        break
      }
    }
    
    if (!is.null(log_det)) {
      if (log_det == -Inf) warning("Laplacian is numerically singular!")
    } else {
      warning(
        "Numerical issues with the determinant. Setting it to zero. Consider increasing `eps`."
      )
      
      log_det <- -Inf
    }
    
    log_ret <- log_ret + log_det
  }
  
  return (log_ret)
}


# Avoid `CholWishart::lmvgamma` from returning an array
lmvgamma <- function (x, p) return (as.double(CholWishart::lmvgamma(x, p)))


# Compute the similarity function based on matrix t-distribution
compute_log_simil_f_t <- function (X, n, p) {
  delta <- 3  # Degrees of freedom of the Hyper Inverse Wishart prior
  delta_star <- 3 + n
  D <- diag(p)
  D_star <- D + crossprod(X)
  
  exp_0 <- 0.5 * (delta + p - 1)
  exp_star <- 0.5 * (delta_star + p - 1)
  
  return (
    lmvgamma(x = exp_star, p = p) - 0.5 * n * p * log(pi) -
      lmvgamma(x = exp_0, p = p) -
      exp_star * as.double(determinant(D_star)$modulus)
  )
}


centers2tessellation <- function (X_cor, centers) {
  # Function that deterministically assigns nodes to centers to obtain a
  # tessellation consisting of supernodes
  # `X_cor` is the matrix with pairwise absolute correlations.
  # `centers` is a vector for, or an indicator vector of variables that are
  # centers.
  p <- NCOL(X_cor)
  assignment <- integer(p)
  tmp <- X_cor[, centers, drop = FALSE]
  for (i in 1:p) assignment[i] <- which.max(tmp[i, ])
  return (assignment)
}


## MCMC update for a graph

# MCMC step based on the Laplace approximation of the normalising constant of
# the G-Wishart distribution
withr::with_makevars(
  new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
  code = Rcpp::sourceCpp(file = "src/laplace.cpp")
)


log_p_GGM_Laplace <- function (adj, df_0, U, n) {
  # Approximate marginal loglikelihood of the principal component scores when
  # the graph is `adj`.
  # `U` is the scatter matrix of the principal component scores.
  adj <- as.matrix(adj)
  p <- nrow(adj)
  rate_0 <- diag(p)
  
  return (log_gwish_norm_laplace_Rcpp(
    adj = adj, df = df_0 + n, rate = rate_0 + U
  ) - 0.5 * n * p * log(2 * pi) - log_gwish_norm_laplace_Rcpp(
    adj = adj, df = df_0, rate = rate_0
  ))
}


update_G_Laplace <- function (adj, edge_prob, U, n, memoise_Laplace) {
  # Update the supergraph using the Laplace approximation for the likelihood.
  # `U` is the scatter matrix of the principal component scores.
  p <- NROW(adj)
  if (p == 1L) return (adj)  # Only one graph possible on one node.
  
  # Sample which edge to flip
  ind <- combn(x = p, m = 2L)[, sample.int(n = choose(p, 2L), size = 1L)]
  adj_prop <- adj
  adj_prop[ind[1], ind[2]] <- 1L - adj[ind[1], ind[2]]
  adj_prop[ind[2], ind[1]] <- adj_prop[ind[1], ind[2]]
  log_R <- (memoise_Laplace(adj_prop, U) - memoise_Laplace(adj, U))
  
  if (edge_prob != 0.5) log_R <- log_R + (2L * adj[ind[1], ind[2]] - 1L) * (
    log1p(-edge_prob) - log(edge_prob)
  )
  
  if (runif(1L) < exp(log_R)) adj <- adj_prop
  return (adj)
}


log_gwish_norm_laplace_empty_graph <- function (p, df, rate_diag) {
  # Log of the "diagonal-Laplace" approximation to the normalising constant of
  # a G-Wishart distribution with an empty graph, degrees of of freedom `df` and
  # rate matrix with diagonal `rate_diag`.
  # Any (additive) terms that do not vary with `df` or `rate_diag` are not
  # included.
  return (0.5 * (p * (-df + (df - 1) * log(df - 2)) - df * sum(log(rate_diag))))
}


log_p_GGM_Laplace_empty_graph <- function (p, df_0, U_diag, n) {
  # Approximate marginal loglikelihood of the principal component scores when
  # the graph is empty.
  # `U_diag` is the diagonal of the scatter matrix of the principal component
  # scores.
  return (log_gwish_norm_laplace_empty_graph(
    p = p, df = df_0 + n, rate_diag = 1 + U_diag
  ) - 0.5 * n * p * log(2 * pi) - log_gwish_norm_laplace_empty_graph(
    p = p, df = df_0, rate_diag = rep(1, p)
  ))
}


# Laplace approximation separately for the unaugmented part and the augmentation
compute_log_p_Y <- function (
  U_unaugmented, U_augmented_diag, G, assignment, n, memoise_Laplace, coars_pow
) {
  # `U` is the scatter matrix of the principal component scores Y.
  # `U_unaugmented` corresponds to the unaugmented `Y_star`.
  # `U_augmented_diag` is the diagonal of the scatter matrix corresponding to
  # the augmentation.
  p <- length(assignment)
  n_centers <- length(unique(assignment))
  G_unaugmented <- as.matrix(G)
  log_p_Y_unaugmented <- memoise_Laplace(adj = G_unaugmented, U = U_unaugmented)
  p_aux <- p - n_centers
  
  log_p_Y_augmentation <- if (p_aux == 0) 0 else log_p_GGM_Laplace_empty_graph(
    p = p_aux, df_0 = df_0, U_diag = U_augmented_diag, n = n
  ) * coars_pow
  
  # Laplace approximation to the marginal likelihood
  return (list(
    log_ret = log_p_Y_unaugmented + log_p_Y_augmentation,
    log_p_Y_unaugmented = log_p_Y_unaugmented
  ))
}


compute_U_unaugmented <- function (X, assignment) {
  # `U_unaugmented` is the scatter matrix of the unaugmented `Y_star`.
  n <- nrow(X)
  n_centers <- max(assignment)
  Y_star <- matrix(data = NA_real_, nrow = n, ncol = n_centers)
  
  for (center in 1:n_centers) {
    ind <- which(assignment == center)
    x_k <- X[, ind, drop = FALSE]
    svd_res <- fast_svd(x = x_k, nu = 0L, n = n, p = length(ind))
    Y_star[, center] <- x_k %*% (svd_res$vt[1, ] / svd_res$d[1])
  }
  
  return (n * crossprod(Y_star))
}


compute_log_posterior_no_G <- function (
    X, X_cor, centers, memoise_tree_activation_f, log_cohesion, nested
) {
  # Part of the posterior that does not vary with the supergraph G to make
  # better use of memoisation
  # `centers` is an indicator vector indicating which nodes are centers.
  n <- NROW(X)
  p <- NCOL(X)
  if (n < p) stop("`n ≥ p` is required.")
  assignment <- centers2tessellation(X_cor, centers)
  p_k <- table(assignment)
  n_centers <- length(p_k)
  log_p_C <- -lchoose(p, n_centers) + sum(log_cohesion[p_k])
  
  for (center in 1:n_centers) {
    ind <- which(assignment == center)
    log_p_C <- log_p_C + memoise_tree_activation_f(ind)
  }
  
  if (is.na(log_p_C)) stop("NA in the computation of the similarities.")
  
  if (nested) {
    U_unaugmented <- NULL
    U_augmented_diag <- NULL
  } else {
    U_unaugmented <- compute_U_unaugmented(X, assignment)
    U_augmented_diag <- n * rep(1, p - n_centers)
  }
  
  return (list(
    log_ret = log_p_C, U_unaugmented = U_unaugmented,
    U_augmented_diag = U_augmented_diag, assignment = assignment, n = n
  ))
}


compute_log_posterior <- function (
    centers, G, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow_lik, nested, pi_se
) {
  # `centers` is an indicator vector indicating which nodes are centers.
  res <- memoise_posterior_no_G(centers, memoise_tree_activation_f)
  
  if (!nested) {
    tmp <- compute_log_p_Y(
      res$U_unaugmented, res$U_augmented_diag, G, res$assignment, res$n,
      memoise_Laplace, coars_pow_lik
    )
    
    res$log_ret <- res$log_ret + tmp$log_ret
    res$log_p_Y <- tmp$log_ret
    res$log_p_Y_unaugmented <- tmp$log_p_Y_unaugmented
    
    # Add prior on supergraph.
    n_centers <- sum(centers)
    n_edges <- sum(G[upper.tri(G)])
    
    res$log_ret <- res$log_ret + n_edges * log(pi_se) +
      (choose(n_centers, 2L) - n_edges) * log(1 - pi_se)
  }
  
  if (is.na(res$log_ret)) stop("log_post is NA!")
  res$assignment <- NULL
  res$U <- res$U_unaugmented
  res$U_unaugmented <- NULL
  return (res)
}


update_G <- function (G, n, tmp_cur, memoise_Laplace, pi_se) {
  G_old <- G
  
  G <- update_G_Laplace(
    adj = G_old, edge_prob = pi_se, U = tmp_cur$U, n, memoise_Laplace
  )
  
  # Update `tmp_cur` accordingly.
  log_p_Y_unaugmented <- memoise_Laplace(G, tmp_cur$U)
  log_diff <- log_p_Y_unaugmented - tmp_cur$log_p_Y_unaugmented
  tmp_cur$log_p_Y_unaugmented <- log_p_Y_unaugmented
  tmp_cur$log_ret <- tmp_cur$log_ret + log_diff
  n_edges <- sum(G[upper.tri(G)])
  n_edges_old <- sum(G_old[upper.tri(G_old)])
  
  # Add term corresponding to the supergraph prior.
  if (n_edges != n_edges_old) tmp_cur$log_ret <-
    tmp_cur$log_ret + (n_edges - n_edges_old) * (log(pi_se) - log(1 - pi_se))
  
  return (list(G = G, tmp_cur = tmp_cur))
}


# Birth-death MCMC step
# This algorithm is based around the addition and deletion of centers to
# generate birth-death steps for supernodes.
take_birth_death_step <- function (
  X_cor, assignment, pi_se, pi_se_prop, centers, G, tmp_cur,
  memoise_tree_activation_f, memoise_posterior_no_G, memoise_Laplace,
  coars_pow_prior, coars_pow_lik, nested
) {
  p <- length(centers)
  n_centers <- sum(centers)
  
  # Sample whether to do a birth or a death
  prob_birth <- 0.5
  if (n_centers == 1L) prob_birth <- 1
  if (n_centers == p) prob_birth <- 0
  birth <- (runif(1L) <= prob_birth)
  
  centers_prop <- centers
  
  if (birth) {  # Birth proposal
    tmp <- sample.int(n = p - n_centers, size = 1L)
    centers_prop[!centers][tmp] <- TRUE
    n_centers_prop <- n_centers + 1
    
    # New assignments
    assignment_prop <- centers2tessellation(X_cor, centers_prop)
    # find index of the new supernode
    sampled_supernode <- sum(centers_prop[1:which(centers_prop != centers)])
    
    # Find which supernodes contain the same variables as before and which have
    # changed.
    supernode_ischanged <- rep(FALSE, n_centers)
    out_h <- outer(assignment_prop, assignment_prop, "==")
    
    for (h in 1:n_centers) {
      ind_h <- which(assignment == h)
      
      if (
        (sum(out_h[ind_h, ind_h] == 0) > 0) ||
          (sum(out_h[ind_h, -ind_h] == 1) > 0)
      ) supernode_ischanged[h] <- TRUE
    }
    
    # Extra supernode goes on the outer edges of adj matrix.
    G_prop <- matrix(NA, nrow = n_centers_prop, ncol = n_centers_prop)
    diag(G_prop) <- 0
    
    # Nodes that did not change
    unchanged_nodes_G <- c(1:n_centers)[!supernode_ischanged]
    
    unchanged_nodes_prop <-
      unchanged_nodes_G + as.numeric(unchanged_nodes_G > sampled_supernode)
    
    if (length(unchanged_nodes_G) > 1L) {
      G_prop[unchanged_nodes_prop,unchanged_nodes_prop] <-
        G[unchanged_nodes_G,unchanged_nodes_G]
    }
    
    # Nodes that changed
    changed_nodes_prop <- setdiff(c(1:n_centers_prop), unchanged_nodes_prop)
    ncount_G_prop <- sum(is.na(G_prop)) / 2
    ncount_G <- ncount_G_prop - n_centers
    nedges_G_prop <- 0
    nedges_G <- 0
    
    if (length(changed_nodes_prop) > 0) {
      for (h1 in changed_nodes_prop) {
        for (h2 in (1:n_centers_prop)[-h1]) {
          if (is.na(G_prop[h1,h2])) {
            G_prop[h1,h2] <- as.numeric(runif(1L) <= pi_se_prop)
            G_prop[h2,h1] <- G_prop[h1,h2]
            nedges_G_prop <- nedges_G_prop + G_prop[h1,h2]
            
            if ((h1 != sampled_supernode) & (h2 != sampled_supernode)) {
              nedges_G <- nedges_G + G[
                h1 - as.numeric(h1 > sampled_supernode),
                h2 - as.numeric(h2 > sampled_supernode)
              ]
            }
          }
        }
      }
    }
    
    # Add discrete uniform for choice of center to add/remove.
    R <- log(p - n_centers) - log(n_centers_prop)
    
    if (!nested) {
      # Proposal term considers the difference in number of edges times logit.
      R <- R + (nedges_G - nedges_G_prop) * log(pi_se_prop) -
        (n_centers + nedges_G - nedges_G_prop) * log(1 - pi_se_prop)
    }
  } else {  # Death proposal
    sampled_supernode <- sample.int(n = n_centers, size = 1L)
    centers_prop[centers][sampled_supernode] <- FALSE
    n_centers_prop <- n_centers - 1L
    
    # New assignments
    assignment_prop <- centers2tessellation(X_cor, centers_prop)
    
    # Find which supernodes contain the same variables as before and which have
    # changed.
    supernode_ischanged <- rep(FALSE, n_centers)
    out_h <- outer(assignment_prop, assignment_prop, "==")
    
    for (h in 1:n_centers) {
      ind_h <- which(assignment == h)
      
      if (
        (sum(out_h[ind_h, ind_h] == 0) > 0) ||
          (sum(out_h[ind_h, -ind_h] == 1) > 0)
      ) supernode_ischanged[h] <- TRUE
    }
    
    # Extra supernode goes on the outer edges of adj matrix.
    G_prop <- matrix(NA, nrow = n_centers_prop, ncol = n_centers_prop)
    diag(G_prop) <- 0
    
    # Nodes that did not change
    unchanged_nodes_G <- c(1:n_centers)[!supernode_ischanged]
    
    unchanged_nodes_prop <-
      unchanged_nodes_G - as.numeric(unchanged_nodes_G > sampled_supernode)
    
    if (length(unchanged_nodes_G) > 1) {
      G_prop[unchanged_nodes_prop, unchanged_nodes_prop] <-
        G[unchanged_nodes_G, unchanged_nodes_G]
    }
    
    # Nodes that changed
    changed_nodes_prop <- setdiff(
      c(1:n_centers)[supernode_ischanged], sampled_supernode
    )
    
    changed_nodes_prop <- changed_nodes_prop - as.numeric(
      changed_nodes_prop > sampled_supernode
    )
    
    ncount_G_prop <- sum(is.na(G_prop)) / 2
    ncount_G <- ncount_G_prop + n_centers_prop
    nedges_G_prop <- 0
    nedges_G <- 0
    
    if (length(changed_nodes_prop) > 0) {
      for (h1 in changed_nodes_prop) {
        for (h2 in (1:n_centers_prop)[-h1]) {
          if (is.na(G_prop[h1,h2])) {
            G_prop[h1,h2] <- as.numeric(runif(1L) <= pi_se_prop)
            G_prop[h2,h1] <- G_prop[h1,h2]
            nedges_G_prop <- nedges_G_prop + G_prop[h1,h2]
            
            nedges_G <- nedges_G + G[
              h1 + as.numeric(h1 >= sampled_supernode),
              h2 + as.numeric(h2 >= sampled_supernode)
            ]
          }
        }
      }
      
      # Count removed part.
      h1 <- sampled_supernode
      for (h2 in (1:n_centers)[-h1]) nedges_G <- nedges_G + G[h1, h2]
    }
    
    # Add discrete uniform for choice of center to add/remove.
    R <- log(n_centers) - log(p - n_centers_prop)
    
    if (!nested) {
      # Proposal term considers the difference in number of edges times logit
      R <- R + (nedges_G - nedges_G_prop) * log(pi_se_prop) -
        (1L - n_centers + nedges_G - nedges_G_prop) * log(1 - pi_se_prop)
    }
    
  }
  
  # Add birth/death probabilities
  if (prob_birth == 1 || prob_birth == 0) {
    R <- R - log(2)
  } else {
    R <- R + (1 - 2 * as.numeric(birth)) * (
      log(prob_birth) - log(1 - prob_birth)
    )
  }
  
  tmp_prop <- compute_log_posterior(
    centers_prop, G_prop, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow_lik, nested, pi_se
  )
  
  if (runif(1L) < exp(R + tmp_prop$log_ret - tmp_cur$log_ret)) {
    assignment <- assignment_prop
    centers <- centers_prop
    G <- G_prop
    tmp_cur <- tmp_prop
  }
  
  return (list(assignment = assignment, centers = centers, G = G, tmp_cur = tmp_cur))
}


# Move MCMC step
take_move_step <- function (
  X_cor, assignment, pi_se, pi_se_prop, centers, G, tmp_cur,
  memoise_tree_activation_f, memoise_posterior_no_G, memoise_Laplace,
  coars_pow_prior, coars_pow_lik, nested
) {
  p <- length(centers)
  n_centers <- sum(centers)
  
  if (n_centers == p) return (
    list(assignment = assignment, centers = centers, G = G, tmp_cur = tmp_cur)
  )
  
  # Center index of node i that is no longer a center in the proposal
  i_c <- sample.int(n = n_centers, size = 1L)
  
  centers_prop <- centers
  centers_prop[centers][i_c] <- FALSE
  
  centers_prop[!centers_prop][sample.int(n = sum(!centers_prop), size = 1L)] <-
    TRUE
  
  # new assignments
  assignment_prop <- centers2tessellation(X_cor, centers_prop)
  
  ass_tab <- table(assignment, assignment_prop)
  
  # Find which supernodes contain the same variables as before and which have
  # changed.
  
  supernode_ischanged <- rep(FALSE, n_centers)
  
  for (h in 1:n_centers) {
    if (sum(ass_tab[h, ] > 0) > 1 || sum(ass_tab[, h] > 0) > 1) {
      supernode_ischanged[h] <- TRUE
    }
  }
  
  # Propose new edges uniformly at random for those supernodes that are
  # different.
  G_prop <- G * NA
  diag(G_prop) <- 0
  
  unchanged_nodes <- c(1:n_centers)[!supernode_ischanged]
  
  if (length(unchanged_nodes) > 1) G_prop[unchanged_nodes,unchanged_nodes] <-
    G[unchanged_nodes,unchanged_nodes]
  
  changed_nodes <- c(1:n_centers)[supernode_ischanged]
  n_changed_nodes <- length(changed_nodes)
  nedges_G_prop <- 0
  nedges_G <- 0
  
  if (n_changed_nodes > 0) {
    for (h1 in changed_nodes) {
      for (h2 in (1:n_centers)[-h1]) {
        if (is.na(G_prop[h1,h2])) {
          G_prop[h1,h2] <- as.numeric(runif(1L) <= pi_se_prop)
          G_prop[h2,h1] <- G_prop[h1,h2]
          nedges_G_prop <- nedges_G_prop + G_prop[h1,h2]
          nedges_G <- nedges_G + G[h1,h2]
        }
      }
    }
  }
  
  tmp_prop <- compute_log_posterior(
    centers_prop, G_prop, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow_lik, nested, pi_se
  )
  
  R <- 0
  
  if (!nested) {
    # Proposal term considers the difference in number of edges times logit
    R <- R + (nedges_G - nedges_G_prop) * (
      log(pi_se_prop) - log(1 - pi_se_prop)
    )
  }
  
  if (runif(1L) < exp(R + tmp_prop$log_ret - tmp_cur$log_ret)) {
    assignment <- assignment_prop
    centers <- centers_prop
    G <- G_prop
    tmp_cur <- tmp_prop
  }
  
  return (
    list(assignment = assignment, centers = centers, G = G, tmp_cur = tmp_cur)
  )
}


spectralClustering <- function (corr_mat, k = 3L) {
  # `k` is the number of clusters.
  # Compute the Laplacian matrix.
  laplacian <- diag(rowSums(corr_mat)) - corr_mat
  
  # Compute the first k eigenvectors of the Laplacian matrix
  eig <- eigen(laplacian, symmetric = TRUE, only.values = FALSE)
  eig_vec <- eig$vectors[, 1:k]
  
  # Normalize the eigenvectors
  norm_eig_vec <- t(apply(eig_vec, 1, function (x) x/sqrt(sum(x^2))))
  
  # Perform k-means clustering on the normalized eigenvectors
  
  # Since kmeans uses a random initialization and we want reproducibility,
  # we set a seed.
  old <- .Random.seed
  set.seed(1L)
  res <- kmeans(norm_eig_vec, k)$cluster
  .Random.seed <<- old
  return (res)
}


get_superedge_mat <- function (assignment, G) {
  p <- length(assignment)
  n_centers <- NROW(G)
  superedge_mat <- matrix(0L, nrow = p, ncol = p)
  if (n_centers == 1L) return (superedge_mat)
  
  for (k in 1:(n_centers - 1L)) {
    ind_k <- which(assignment == k)
    
    for (l in (k + 1L):n_centers) {
      ind_l <- which(assignment == l)
      superedge_mat[ind_k, ind_l] <- G[k, l]
      superedge_mat[ind_l, ind_k] <- G[k, l]
    }
  }
  
  return (superedge_mat)
}


# Main function that runs the MCMC
# `nested` indicates whether nested MCMC or the algorithm with coarsening of the
# likelihood should be preformed.
# `n_iter` is the number of MCMC iterations after burn-in.
# `log_cohesion` is the log of the cohesion function.
# `burnin` is the number of iterations discarded as burn-in,
# which equals `n_iter %/% 10L` by default.
# `coars_pow` is the coarsening parameter ($\zeta$ in the text).
# `pi_se` is the prior superedge inclusion probability.
# `pi_se_prop` is the superedge probability in the Metropolis-Hasting proposal.
# `use_t` indicates whether to use a similarity function based on the matrix t
# distribution.
run_MCMC <- function (
  X, nested, n_iter, burnin = NULL, n_thin = n_thin, log_cohesion,
  coars_pow = NULL, pi_se = NULL, pi_se_prop = NULL, use_t = FALSE
) {
  p <- NCOL(X)
  n <- NROW(X)
  if (is.null(burnin)) burnin <- n_iter %/% 10L
  if (is.null(coars_pow)) coars_pow <- 10 / n
  if (is.null(pi_se)) pi_se <- 0.5
  if (is.null(pi_se_prop)) pi_se_prop <- 0.5
  
  coars_pow_prior <- coars_pow
  coars_pow_lik <- coars_pow
  X_cor <- abs(cor(X))
  log_w_matrix <- compute_log_w_matrix(X, coars_pow_prior)
  centers_MCMC <- matrix(data = NA, nrow = n_iter, ncol = p)
  superedge_mat <- matrix(0L, nrow = p, ncol = p)
  post_sim <- matrix(data = 0L, nrow = p, ncol = p)
  log_post_MCMC <- numeric(n_iter)
  G_MCMC <- vector(mode = "list", length = n_iter)
  
  memoise_tree_activation_f <- memoise::memoise(f = if (use_t) function (ind) {
    coars_pow * compute_log_simil_f_t(
      X = X[, ind, drop = FALSE], n, length(ind)
    )
  } else function (ind) {
    compute_log_tree_activation_f(
      log_w_matrix = log_w_matrix[ind, ind], p = length(ind)
    )
  }, cache = cachem::cache_mem(max_size = Inf))
  
  memoise_posterior_no_G <- memoise::memoise(
    f = function (centers, memoise_tree_activation_f) {
      compute_log_posterior_no_G(
        X, X_cor, centers, memoise_tree_activation_f, log_cohesion, nested
      )
    }, cache = cachem::cache_mem(max_size = Inf),
    omit_args = "memoise_tree_activation_f"
  )
  
  memoise_Laplace <- memoise::memoise(f = function (adj, U) {
    coars_pow_lik * log_p_GGM_Laplace(adj, df_0, U, n)
  }, cache = cachem::cache_mem(max_size = Inf))
  
  K_init <- 3L
  centers <- logical(p)
  
  if (p < K_init) centers[] <- TRUE else {
    # Initialize using spectral clustering.
    centers[match(1:K_init,  spectralClustering(X_cor, K_init))] <- TRUE
  }
  
  n_centers <- sum(centers)
  
  # Initilise assignments
  assignment <- centers2tessellation(X_cor, centers)
  
  
  # Supergraph: initialize at empty
  G <- matrix(0L, nrow = n_centers, ncol = n_centers)
  
  tmp_cur <- compute_log_posterior(
    centers, G, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow_lik, nested, pi_se
  )
  
  pb <- txtProgressBar(max = burnin + n_iter, style = 3)
  
  for (s in 1:(burnin + n_iter)) {
    tmp <- take_birth_death_step(
      X_cor, assignment, pi_se, pi_se_prop, centers, G, tmp_cur,
      memoise_tree_activation_f, memoise_posterior_no_G, memoise_Laplace,
      coars_pow_prior, coars_pow_lik, nested
    )
    
    assignment <- tmp$assignment
    tmp_cur <- tmp$tmp_cur
    centers <- tmp$centers
    G <- tmp$G
    
    if (!nested) {
      # Update G by itself.
      res <- update_G(G, n, tmp_cur, memoise_Laplace, pi_se)
      G <- res$G
      tmp_cur <- res$tmp_cur
    }
    
    tmp <- take_move_step(
      X_cor, assignment, pi_se, pi_se_prop, centers, G, tmp_cur,
      memoise_tree_activation_f, memoise_posterior_no_G, memoise_Laplace,
      coars_pow_prior, coars_pow_lik, nested
    )
    
    assignment <- tmp$assignment
    tmp_cur <- tmp$tmp_cur
    centers <- tmp$centers
    G <- tmp$G
    
    if (!nested) {
      # Update G by itself.
      res <- update_G(G, n, tmp_cur, memoise_Laplace, pi_se)
      G <- res$G
      tmp_cur <- res$tmp_cur
    }
    
    if (s > burnin) {
      centers_MCMC[s - burnin, ] <- tmp$centers
      log_post_MCMC[s - burnin] <- tmp_cur$log_ret
      G_MCMC[[s - burnin]] <- G
      post_sim <- post_sim + outer(assignment, assignment, "==")
      
      if (!nested) {
        # Keep track of the superedges.
        superedge_mat <- superedge_mat + get_superedge_mat(assignment, G)
      }
    }
    
    setTxtProgressBar(pb, s)
  }
  
  close(pb)
  
  if (nested) {
    # Inner MCMC
    # Compute the superedge inclusion probabilities for the nested MCMC.
    n_recorded <- n_iter %/% n_thin
    superedge_mat <- matrix(0, nrow = p, ncol = p)
    print("Running inner MCMC...")
    set.seed(1L)
    pb <- txtProgressBar(max = n_recorded, style = 3)
    
    for (s in 1:n_recorded) {
      centers <- centers_MCMC[s * n_thin, ]
      n_centers <- sum(centers)
      assignment <- centers2tessellation(X_cor, centers)
      U <- compute_U_unaugmented(X, assignment)
      
      # Supergraph: initialize to empty graph.
      G <- matrix(0L, nrow = n_centers, ncol = n_centers)
      
      for (s_inner in 1:1000) G <- update_G_WWA(G, pi_se, 3, U, n, 1)
      superedge_mat <- superedge_mat + get_superedge_mat(assignment, G)
      setTxtProgressBar(pb, s)
    }
    
    close(pb)
  }
  
  return (list(
    centers_MCMC = centers_MCMC,
    superedge_mat = superedge_mat / if (nested) n_recorded else n_iter,
    post_sim = post_sim / n_iter,
    log_post_MCMC = log_post_MCMC,
    G_MCMC = G_MCMC
  ))
}
