### Graph sphere: MCMC

# Implementation of the Markov chain Monte Carlo described in the paper "Graph
# Sphere: From Nodes to Supernodes in Graphical Models" by Willem van den Boom,
# Maria De Iorio, Alexandros Beskos and Ajay Jasra.



# Install required packages
pkgs <- c(
  "abind", "BH", "latticeExtra", "matrixStats", "memoise", "mvtnorm", "qgam",
  "Rcpp", "R.utils", "withr"
)

pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(pkgs) > 0L) {
  print("Installing required packages...")
  
  install.packages(
    pkgs = pkgs, dependencies = TRUE, repos = "https://cloud.r-project.org"
  )
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
update_G_WWA <- function (adj, edge_prob, df_0, U, n) {
  p <- nrow(adj)
  
  return (tryCatch(
    expr = update_G_Rcpp(
      adj = adj,
      
      edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
      df = df_0 + n,
      df_0 = df_0,
      rate = diag(p) + U,
      n_edge = p,  # Number of single edge updates attempted in this MCMC update
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
  if(is.null(n)) n <- nrow(x)
  if(is.null(p)) p <- ncol(x)
  
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
  # Any terms that can be absorbed in the normalising constant of the
  # tessellation prior are dropped. Thus, the function computes the log of
  # $p_k^{2 - p_k} \sum_{E_k} \prod_{(i, j)\in E_k} w_{ij}$.
  log_diag_Laplacian <- numeric(p)
  
  for (i in 1:p) log_diag_Laplacian[i] <- if (p == 1L) -Inf else {
    matrixStats::logSumExp(log_w_matrix[i, -i])
  }
  
  log_ret <- (2L - p) * log(p) + sum(diag(as.matrix(log_w_matrix)))
  if (p == 2L) log_ret <- log_ret + log_diag_Laplacian[1]
  
  if (p > 2L) {
    
    # Cycle through all Laplacian minors to maximize numerical stability.
    # Only if !use_fast_svd
    n_count <- 1L
    log_det <- NULL
    count_det <- 0L
    Laplacian <- diag(p)
    upper_ind <- upper.tri(Laplacian)
    
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
  # better use of memoization
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
  
  # Add p(G | C) term, uniform over all possible supergraphs
  log_ret <- log_p_C - if (nested) 0 else {
    ((n_centers * (n_centers - 1L)) %/% 2L) * log(2)
  }
  
  return (list(
    log_ret = log_ret, U_unaugmented = U_unaugmented,
    U_augmented_diag = U_augmented_diag, assignment = assignment, n = n
  ))
}


compute_log_posterior <- function (
  centers, G, memoise_tree_activation_f, memoise_posterior_no_G,
  memoise_Laplace, coars_pow, nested
) {
  # `centers` is an indicator vector indicating which nodes are centers.
  res <- memoise_posterior_no_G(centers, memoise_tree_activation_f)
  
  if (!nested) {
    tmp <- compute_log_p_Y(
      res$U_unaugmented, res$U_augmented_diag, G, res$assignment, res$n,
      memoise_Laplace, coars_pow
    )
    
    res$log_ret <- res$log_ret + tmp$log_ret
    res$log_p_Y <- tmp$log_ret
    res$log_p_Y_unaugmented <- tmp$log_p_Y_unaugmented
  }
  
  if (is.na(res$log_ret)) stop("log_post is NA!")
  res$assignment <- NULL
  res$U <- res$U_unaugmented
  res$U_unaugmented <- NULL
  return (res)
}


update_G <- function (G, n, tmp_cur, memoise_Laplace) {
  G <- update_G_Laplace(
    adj = G, edge_prob = 0.5, U = tmp_cur$U, n, memoise_Laplace
  )
  
  # Update `tmp_cur` accordingly.
  log_p_Y_unaugmented <- memoise_Laplace(G, tmp_cur$U)
  log_diff <- log_p_Y_unaugmented - tmp_cur$log_p_Y_unaugmented
  tmp_cur$log_p_Y_unaugmented <- log_p_Y_unaugmented
  tmp_cur$log_ret <- tmp_cur$log_ret + log_diff
  return (list(G = G, tmp_cur = tmp_cur))
}


# Birth-death MCMC step
# This algorithm is based around the addition and deletion of centers to
# generate birth-death steps for supernodes.
take_birth_death_step <- function (
  centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
  memoise_Laplace, coars_pow, nested
) {
  p <- length(centers)
  n_centers <- sum(centers)
  
  # Sample whether to do a birth or a death
  prob_birth <- 0.5
  if (n_centers == 1L) prob_birth <- 1
  if (n_centers == p) prob_birth <- 0
  birth <- runif(1L) < prob_birth
  centers_prop <- centers
  
  if (birth) {  # Birth proposal
    tmp <- sample.int(n = p - n_centers, size = 1L)
    centers_prop[!centers][tmp] <- TRUE
    sampled_supernode <- sum(centers_prop[1:which(centers_prop != centers)])
    G_prop <- matrix(0L, nrow = n_centers + 1L, ncol = n_centers + 1L)
    G_prop[-sampled_supernode, -sampled_supernode] <- G
    
    G_prop[sampled_supernode, -sampled_supernode] <- rbinom(
      n = n_centers, size = 1L, prob = 0.5
    )
    
    G_prop[-sampled_supernode, sampled_supernode] <- G_prop[
      sampled_supernode, -sampled_supernode
    ]
    
    R <- (ifelse(
      n_centers + 1L == p, 1, 1 / 2
    ) / (n_centers + 1L)) / (prob_birth / (p - n_centers))
    
    # Add proposal term related to G.
    if (!nested) R <- R / 0.5^n_centers
  } else {  # Death proposal
    sampled_supernode <- sample.int(n = n_centers, size = 1L)
    centers_prop[centers][sampled_supernode] <- FALSE
    G_prop <- G[-sampled_supernode, -sampled_supernode, drop = FALSE]
    
    R <- (ifelse(
      n_centers == 2L, 1, 1 / 2
    ) / (p - n_centers + 1L)) / ((1 - prob_birth) / n_centers)
    
    # Add proposal term related to G.
    if (!nested) R <- R * 0.5^(n_centers - 1L)
  }
  
  tmp_prop <- compute_log_posterior(
    centers_prop, G_prop, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow, nested
  )
  
  if(runif(1) < R * exp(tmp_prop$log_ret - tmp_cur$log_ret)) {
    centers <- centers_prop
    G <- G_prop
    tmp_cur <- tmp_prop
  }
  
  return (list(centers = centers, G = G, tmp_cur = tmp_cur))
}


# Move MCMC step
take_move_step <- function (
    centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow, nested
) {
  p <- length(centers)
  n_centers <- sum(centers)
  
  if (n_centers == p) {
    return (list(centers = centers, G = G, tmp_cur = tmp_cur))
  }
  
  # Center index of node i that is no longer a center in the proposal
  i_c <- sample.int(n = n_centers, size = 1L)
  
  centers_prop <- centers
  centers_prop[centers][i_c] <- FALSE
  centers_prop[!centers][sample.int(n = p - n_centers, size = 1L)] <- TRUE
  
  # Center index of node j that is a new center in the proposal
  j_c <- sum(centers_prop[1:which(centers_prop & !centers)])
  
  # Relabel the supernodes in the supergraph according to the proposed change to
  # `centers`.
  relabel <- if (i_c == j_c) 1:n_centers else if (i_c < j_c) c(
    seq_len(i_c - 1L), j_c, i_c + seq_len(j_c - 1L - i_c), i_c,
    j_c + seq_len(n_centers - j_c)
  ) else c(
    seq_len(j_c - 1L), i_c, j_c + seq_len(i_c - 1L - j_c), j_c,
    i_c + seq_len(n_centers - i_c)
  )
  
  G_prop <- G[relabel, relabel, drop = FALSE]
  
  tmp_prop <- compute_log_posterior(
    centers_prop, G_prop, memoise_tree_activation_f,
    memoise_posterior_no_G, memoise_Laplace, coars_pow, nested
  )
  
  if(runif(1) < exp(tmp_prop$log_ret - tmp_cur$log_ret)) {
    centers <- centers_prop
    G <- G_prop
    tmp_cur <- tmp_prop
  }
  
  return (list(centers = centers, G = G, tmp_cur = tmp_cur))
}


## The following functions are used to implement the locally balanced informed
## proposal.

compute_log_base_proposal <- function (n_centers, n_centers_prop, p, k, move) {
  # If `move = TRUE`, then only consider moves.
  # If `move = FALSE`, then only consider birth/death steps.
  log_prob_bdm <- rep(-log(3), 3L)
  if (n_centers == p) log_prob_bdm <- c(-Inf, 0, -Inf)
  if (isTRUE(move)) log_prob_bdm <- c(-Inf, -Inf, 0)
  
  if (n_centers == 1L) log_prob_bdm <- if (isFALSE(move)) {
    c(0, -Inf, -Inf)
  } else log(c(0.5, 0, 0.5))
  
  if (n_centers_prop == n_centers + 1L) {  # Birth
    log_n_G_birth <- if (k >= n_centers) n_centers * log(2) else {
      matrixStats::logSumExp(lchoose(n_centers, 0:k))
    }
    
    return (log_prob_bdm[1] - log(p - n_centers) - log_n_G_birth)
  }
  
  if (n_centers_prop == n_centers - 1L) {  # Death
    return (log_prob_bdm[2] - log(n_centers))
  }
  
  if (n_centers_prop == n_centers) {  # Move
    return (log_prob_bdm[3] - log(n_centers) - log(p - n_centers))
  }
  
  return (-Inf)
}


ind2prop <- function (prop_ind, centers, G, k, move = NA) {
  # If `move = TRUE`, then only consider moves.
  # If `move = FALSE`, then only consider birth/death steps.
  p <- length(centers)
  n_centers <- sum(centers)
  
  # Birth
  n_G_birth <- if (isTRUE(move) | n_centers == p) 0L else {
    as.integer(if (k == Inf) 2^n_centers else sum(choose(n_centers, 0:k)))
  }
  
  n_birth <- (p - n_centers) * n_G_birth
  
  if (prop_ind <= n_birth) {  # Birth
    ind_add <- (prop_ind - 1L) %/% n_G_birth + 1L
    G_prop_ind <- (prop_ind - 1L) %% n_G_birth + 1L
    centers_prop <- centers
    centers_prop[!centers][ind_add] <- TRUE
    new_supernode <- sum(centers_prop[1:which(centers_prop != centers)])
    G_prop <- matrix(0L, nrow = n_centers + 1L, ncol = n_centers + 1L)
    G_prop[-new_supernode, -new_supernode] <- G
    
     if (k == Inf) G_prop[new_supernode, -new_supernode][min(
      n_centers - ceiling(log2(G_prop_ind) - 1), n_centers
    ):n_centers] <- as.integer(
      strsplit(R.utils::intToBin(G_prop_ind - 1L), split = "")[[1]]
    ) else if (k == 1L) {
      if (G_prop_ind > 1L) {
        G_prop[new_supernode, -new_supernode][G_prop_ind - 1L] <- 1L
      }
    } else stop("Only k = 1,Inf have been implemented.")
    
    G_prop[-new_supernode, new_supernode] <- G_prop[
      new_supernode, -new_supernode
    ]
    
    return (list(centers_prop = centers_prop, G_prop = G_prop))
  }
  
  prop_ind <- prop_ind - n_birth
  n_death <- if (isTRUE(move) | n_centers == 1L) 0L else {
    death_candidates <- which(colSums(G) <= k)
    length(death_candidates)
  }
  
  if (prop_ind <= n_death) {  # Death
    ind_del <- death_candidates[prop_ind]
    centers_prop <- centers
    centers_prop[centers][ind_del] <- FALSE
    G_prop <- G[-ind_del, -ind_del, drop = FALSE]
    return (list(centers_prop = centers_prop, G_prop = G_prop))
  }
  
  prop_ind <- prop_ind - n_death
  n_move <- if (isFALSE(move)) 0L else n_centers * (p - n_centers)
  
  if (prop_ind <= n_move) {  # Move
    ind_del <- (prop_ind - 1L) %/% (p - n_centers) + 1L
    ind_add <- (prop_ind - 1L) %% (p - n_centers) + 1L
    centers_prop <- centers
    centers_prop[centers][ind_del] <- FALSE
    centers_prop[!centers][ind_add] <- TRUE
    return (list(centers_prop = centers_prop, G_prop = G))
  }
  
  stop("`prop_ind` is too large.")
}


prop2ind <- function (centers_prop, G_prop, centers, G, k, move) {
  # If `move = TRUE`, then only consider moves.
  p <- length(centers)
  n_centers <- sum(centers)
  diff_ind <- which(centers_prop != centers)
  
  n_G_birth <- if (isTRUE(move) | n_centers == p) 0L else {
    as.integer(if (k == Inf) 2^n_centers else sum(choose(n_centers, 0:k)))
  }
  
  if (length(diff_ind) == 1L & !centers[diff_ind[1]]) {  # Birth
    ind_add <- which(centers_prop[!centers] != centers[!centers])
    new_supernode <- sum(centers_prop[1:which(centers_prop != centers)])
    
    G_prop_ind <- if (k == Inf) strtoi(
      x = paste(G_prop[-new_supernode, new_supernode], collapse = ""), base = 2L
    ) + 1L else if (k == 1L) {
      tmp <- G_prop[-new_supernode, new_supernode] == 1L
      if (any(tmp)) which(tmp) + 1L else 1L
    } else stop("Only k = 1,Inf implemented so far.")
    
    return (n_G_birth * (ind_add - 1L) + G_prop_ind)
  }
  
  n_birth <- (p - n_centers) * n_G_birth
  
  if (length(diff_ind) == 1L & centers[diff_ind[1]]) {  # Death
    ind_del <- which(!centers_prop[centers])
    
    ind_del <- ind_del - sum(
      colSums(G[, seq_len(ind_del - 1L), drop = FALSE]) > k
    )
    
    return (n_birth + ind_del)
  }
  
  n_death <- if (isTRUE(move) | n_centers == 1L) 0L else {
    n_centers - sum(colSums(G) > k)
  }
  
  if (length(diff_ind) == 2L & n_centers == sum(centers_prop)) {  # Move
    ind_add <- which(centers_prop[!centers] != centers[!centers])
    ind_del <- which(centers_prop[centers] != centers[centers])
    return (n_birth + n_death + (p - n_centers) * (ind_del - 1L) + ind_add)
  }
  
  stop("Proposed values are not in the neighborhood of the current values.")
}


compute_log_prod <- function (
  centers, G, memoise_tree_activation_f, memoise_posterior_no_G,
  memoise_Laplace, coars_pow, nested, k, move
) {
  # Locally balanced informed proposal
  # This function computes the log of the product of the posterior and the ratio
  # of base proposals.
  # Only supernodes with up to `k` superedges can be birthed or killed.
  # If `move = TRUE`, then only consider moves.
  # If `move = FALSE`, then only consider birth/death steps.
  p <- length(centers)
  n_centers <- sum(centers)
  
  # Birth
  log_prop_ratio_birth <- compute_log_base_proposal(
    n_centers + 1L, n_centers, p, k, move
  ) - compute_log_base_proposal(n_centers, n_centers + 1L, p, k, move)
  
  n_G_birth <- if (isTRUE(move) | n_centers == p) 0L else {
    as.integer(if (k == Inf) 2^n_centers else sum(choose(n_centers, 0:k)))
  }
  
  n_birth <- (p - n_centers) * n_G_birth
  log_prod_birth <- numeric(n_birth)  # Product of posterior and base proposal
  
  for (ind in seq_len(n_birth)) {
    tmp <- ind2prop(ind, centers, G, k)
    
    log_prod_birth[ind] <- log_prop_ratio_birth + compute_log_posterior(
      tmp$centers_prop, tmp$G_prop, memoise_tree_activation_f,
      memoise_posterior_no_G, memoise_Laplace, coars_pow, nested
    )$log_ret
  }
  
  # Death
  log_prop_ratio_death <- compute_log_base_proposal(
    n_centers - 1L, n_centers, p, k, move
  ) - compute_log_base_proposal(n_centers, n_centers - 1L, p, k, move)
  
  n_death <- if (isTRUE(move) | n_centers == 1L) 0L else sum(colSums(G) <= k)
  log_prod_death <- numeric(n_death)
  
  for (ind in seq_len(n_death)) {
    tmp <- ind2prop(n_birth + ind, centers, G, k)
    
    log_prod_death[ind] <- log_prop_ratio_death + compute_log_posterior(
      tmp$centers_prop, tmp$G_prop, memoise_tree_activation_f,
      memoise_posterior_no_G, memoise_Laplace, coars_pow, nested
    )$log_ret
  }
  
  # Move
  n_move <- if (isFALSE(move)) 0L else n_centers * (p - n_centers)
  log_prod_move <- numeric(n_move)
  
  for (ind in seq_len(n_move)) {
    tmp <- ind2prop(ind, centers, G, k, move = TRUE)
    
    log_prod_move[ind] <- compute_log_posterior(
      tmp$centers_prop, G, memoise_tree_activation_f, memoise_posterior_no_G,
      memoise_Laplace, coars_pow, nested
    )$log_ret
  }
  
  return (list(
    log_prod = c(log_prod_birth, log_prod_death, log_prod_move),
    n_bdm = c(n_birth, n_death, n_move)
  ))
}


# Log of the balancing function t / (1 + t)
log_bal_fun <- function (log_t) log_t - qgam::log1pexp(log_t)


take_MCMC_step_informed <- function (
  centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
  memoise_Laplace, coars_pow, nested, move
) {
  # Locally balanced informed proposal
  # If `move = TRUE`, then only consider moves.
  # If `move = FALSE`, then only consider birth/death steps.
  p <- length(centers)
  n_centers <- sum(centers)
  
  # Sample from the informed proposal
  if (move & n_centers == p) {
    return (list(centers = centers, G = G, tmp_cur = tmp_cur))
  }
  
  # Maximum number of superedges affected by birth/death moves
  # ($\kappa$ in the text)
  k <- 1L
  
  log_prod_cur <- compute_log_prod(
    centers, G, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow, nested, k, move
  )
  
  log_prop <- log_bal_fun(log_prod_cur$log_prod - tmp_cur$log_ret) + rep(c(
    compute_log_base_proposal(n_centers, n_centers + 1L, p, k, move),  # Birth
    compute_log_base_proposal(n_centers, n_centers - 1L, p, k, move),  # Death
    compute_log_base_proposal(n_centers, n_centers, p, k, move)  # Move
  ), log_prod_cur$n_bdm)
  
  log_prop <- log_prop - matrixStats::logSumExp(log_prop)
  prop_ind <- sample.int(n = length(log_prop), size = 1L, prob = exp(log_prop))
  prop_list <- ind2prop(prop_ind, centers, G, k, move)
  centers_prop <- prop_list$centers_prop
  G_prop <- prop_list$G_prop
  
  # Evaluate the reverse probability
  tmp_prop <- compute_log_posterior(
    centers_prop, G_prop, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow, nested
  )
  
  n_centers_prop <- sum(centers_prop)
  
  log_prod_prop <- compute_log_prod(
    centers_prop, G_prop, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow, nested, k, move
  )
  
  log_prop_rev <- log_bal_fun(
    log_prod_prop$log_prod - tmp_prop$log_ret
  ) + rep(c(
    # Birth
    compute_log_base_proposal(n_centers_prop, n_centers_prop + 1L, p, k, move),
    # Death
    compute_log_base_proposal(n_centers_prop, n_centers_prop - 1L, p, k, move),
    # Move
    compute_log_base_proposal(n_centers_prop, n_centers_prop, p, k, move)
  ), log_prod_prop$n_bdm)
  
  log_prop_rev <- log_prop_rev - matrixStats::logSumExp(log_prop_rev)
  cur_ind <- prop2ind(centers, G, centers_prop, G_prop, k, move)
  
  # Metropolis-Hastings step
  if (runif(1) < exp(
    tmp_prop$log_ret + log_prop_rev[cur_ind] - tmp_cur$log_ret -
      log_prop[prop_ind]
  )) {
    centers <- centers_prop
    G <- G_prop
    tmp_cur <- tmp_prop
  }
  
  return (list(centers = centers, G = G, tmp_cur = tmp_cur))
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
# `prob_informed` is the probability at each MCMC iteration that the locally
# balanced informed proposals are used.
# `burnin` is the number of iterations discarded as burn-in,
# which equals `n_iter %/% 10L` by default.
# `coars_pow` is the coarsening parameter ($\zeta$ in the text).
run_MCMC <- function (
  X, nested, n_iter, log_cohesion, prob_informed = 0, burnin = NULL,
  coars_pow = NULL
) {
  p <- NCOL(X)
  n <- NROW(X)
  if (is.null(burnin)) burnin <- n_iter %/% 10L
  if (is.null(coars_pow)) coars_pow <- 10 / n
  X_cor <- abs(cor(X))
  log_w_matrix <- compute_log_w_matrix(X, coars_pow)
  centers_MCMC <- matrix(data = NA, nrow = n_iter, ncol = p)
  superedge_mat <- matrix(0L, nrow = p, ncol = p)
  post_sim <- matrix(data = 0L, nrow = p, ncol = p)
  log_post_MCMC <- numeric(n_iter)
  G_MCMC <- vector(mode = "list", length = n_iter)
  
  memoise_tree_activation_f <- memoise::memoise(f = function (ind) {
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
    coars_pow * log_p_GGM_Laplace(adj, df_0, U, n)
  }, cache = cachem::cache_mem(max_size = Inf))
  
  K_init <- 3L
  centers <- logical(p)
  
  if (p < K_init) centers[] <- TRUE else {
    # Initialize using spectral clustering.
    centers[match(1:K_init,  spectralClustering(X_cor, K_init))] <- TRUE
  }
  
  n_centers <- sum(centers)
  
  # Supergraph: initialize at empty
  G <- matrix(0L, nrow = n_centers, ncol = n_centers)
  
  tmp_cur <- compute_log_posterior(
    centers, G, memoise_tree_activation_f, memoise_posterior_no_G,
    memoise_Laplace, coars_pow, nested
  )
  
  pb <- txtProgressBar(max = burnin + n_iter, style = 3)
  
  for (s in 1:(burnin + n_iter)) {
    informed <- runif(1) < prob_informed
    
    tmp <- if (informed) {
      # Locally balanced informed proposal
      tmp <- take_MCMC_step_informed(
        centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
        memoise_Laplace, coars_pow, nested, move = FALSE
      )
      
      centers <- tmp$centers
      G <- tmp$G
      tmp_cur <- tmp$tmp_cur
      
      take_MCMC_step_informed(
        centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
        memoise_Laplace, coars_pow, nested, move = TRUE
      )
    } else {
      tmp <- take_birth_death_step(
        centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
        memoise_Laplace, coars_pow, nested
      )
      
      tmp_cur <- tmp$tmp_cur
      centers <- tmp$centers
      G <- tmp$G
      
      take_move_step(
        centers, G, tmp_cur, memoise_tree_activation_f, memoise_posterior_no_G,
        memoise_Laplace, coars_pow, nested
      )
    }
    
    tmp_cur <- tmp$tmp_cur
    centers <- tmp$centers
    G <- tmp$G
    
    if (!nested) {
      # Update G by itself.
      res <- update_G(G, n, tmp_cur, memoise_Laplace)
      G <- res$G
      tmp_cur <- res$tmp_cur
    }
    
    if (s > burnin) {
      centers_MCMC[s - burnin, ] <- tmp$centers
      log_post_MCMC[s - burnin] <- tmp_cur$log_ret
      G_MCMC[[s - burnin]] <- G
      assignment <- centers2tessellation(X_cor, centers)
      
      for (i in 1:p) {
        post_sim[i, ] <- post_sim[i, ] + (assignment[i] == assignment)
      }
      
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
    superedge_prob <- 0.5
    n_recorded <- min(1e3L, n_iter %/% 10L)
    n_thin <- n_iter %/% n_recorded
    superedge_mat <- matrix(0L, nrow = p, ncol = p)
    print("Running inner MCMC...")
    set.seed(1L)
    pb <- txtProgressBar(max = n_recorded, style = 3)
    
    for (s in 1:n_recorded) {
      centers <- centers_MCMC[s * n_thin, ]
      n_centers <- sum(centers)
      if (n_centers == 1L) next
      assignment <- centers2tessellation(X_cor, centers)
      U <- compute_U_unaugmented(X, assignment)
      
      # Supergraph: initialize at prior draw
      G <- matrix(0L, nrow = n_centers, ncol = n_centers)
      
      G[upper.tri(G)] <- rbinom(
        n = (n_centers * (n_centers - 1L)) %/% 2L, size = 1L,
        prob = superedge_prob
      )
      
      G <- G + t(G)
      for (s_inner in 1:1e3L) G <- update_G_WWA(G, superedge_prob, 3, U, n)
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
