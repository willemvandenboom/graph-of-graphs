### Graph of graphs: Simulation study on inferring individual edges

# Code for the simulation study on inferring individual edges in the paper
# "Graph of Graphs: From Nodes to Supernodes in Graphical Models" by Maria De
# Iorio, Willem van den Boom, Alexandros Beskos, Ajay Jasra and Andrea
# Cremaschi.


# Install required packages
pkgs <- c(
  "doRNG", "doSNOW", "GeneNet", "mvtnorm", "RcppBlaze", "viridisLite", "withr"
)

pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(pkgs) > 0L) {
  print("Installing required packages...")
  install.packages(pkgs = pkgs, repos = "https://cloud.r-project.org")
}



# The following uses the Rcpp function `rgwish_Rcpp` to provide a sampler for
# the G-Wishart distribution.
withr::with_makevars(
  new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
  code = Rcpp::sourceCpp(file = "src/ggm.cpp")
)


#' Sample from the G-Wishart distribution
#' 
#' @description Sample from the G-Wishart distribution following the method of
#' Lenkoski (2013, Section 2.4, arXiv:1304.1350v1).
#' 
#' @param adj The symmetric adjacency matrix of the graph with zeroes on the
#'   diagonal.
#' @param df The degrees of freedom of the G-Wishart distribution.
#' @param rate The inverse of the scale matrix of the G-Wishart distribution.
#' @return A positive definite matrix
rgwish <- function(adj, df = 3, rate = NULL) {
  p <- nrow(adj)
  if (is.null(rate)) rate <- diag(p)
  
  return(rgwish_Rcpp(
    adj = adj,
    df = df,
    rate = rate,
    seed = sample.int(n = .Machine$integer.max, size = 1)
  ))
}


# Generate data for p = 20, 40.
set.seed(1L)
n_rep <- 10L
prop_vec <- c(0, 0.25, 0.5, 0.75)  # Proportion of edges in the true graph
n_prop <- length(prop_vec)
p_vec <- c(20L, 40L)
n_p <- length(p_vec)
n_vec <- c(1e1L, 5e1L, 1e2L, 5e2L, 1e3L)
n_n <- length(n_vec)
G_true_list <- vector("list", n_p)
data_list <- vector("list", n_p)

for (p_ind in 1:n_p) {
  p <- p_vec[p_ind]
  
  G_true_list[[p_ind]] <- array(
    data = NA_integer_, dim = c(n_prop, n_rep, p, p)
  )
  
  data_list[[p_ind]] <- array(
    data = NA_real_, dim = c(n_prop, n_rep, n_vec[n_n], p)
  )
  
  for (r in 1:n_rep) for (prop_ind in 1:n_prop) {
    pcor <- GeneNet::ggm.simulate.pcor(num.nodes = p, etaA = prop_vec[prop_ind])
    G_true <- pcor != 0
    diag(G_true) <- FALSE
    mode(G_true) <- "integer"
    G_true_list[[p_ind]][prop_ind, r, , ] <- G_true
    data <- mvtnorm::rmvnorm(n = n_vec[n_n], sigma = solve(rgwish(G_true)))
    data_list[[p_ind]][prop_ind, r, , ] <- data
  }
}


# Generate data for p = 100.
set.seed(1L)
p <- 100L
G_true_arr <- array(data = NA_integer_, dim = c(n_prop, n_rep, p, p))
data_arr <- array(data = NA_real_, dim = c(n_prop, n_rep, n_vec[n_n], p))

for (r in 1:n_rep) for (prop_ind in 1:n_prop) {
  pcor <- GeneNet::ggm.simulate.pcor(num.nodes = p, etaA = prop_vec[prop_ind])
  G_true <- pcor != 0
  diag(G_true) <- FALSE
  mode(G_true) <- "integer"
  G_true_arr[prop_ind, r, , ] <- G_true
  data <- mvtnorm::rmvnorm(n = n_vec[n_n], sigma = solve(rgwish(G_true)))
  data_arr[prop_ind, r, , ] <- data
}


# Remove Rcpp functions as they are not linked properly in the parallel workers
# otherwise.
rm(rgwish_Rcpp)
rm(update_G_Rcpp)


n_iter <- 1e5L
n_cores <- parallel::detectCores()
print(paste(n_cores, "cores detected"))


# Compute results for p = 20, 40 using parallel computing.

run_task <- function (prop_ind, r, n_ind, p_ind) {
  withr::with_makevars(
    new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
    code = Rcpp::sourceCpp(file = "src/ggm.cpp")
  )
  
  
  ## MCMC update for a graph
  
  # We now use the function `update_G_Rcpp` to define an R function that
  # provides an MCMC update for a graph in a Gaussian graphical model.
  
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
  update_G <- function(adj, edge_prob, df_0, U, n) {
    p <- nrow(adj)
    
    return (tryCatch(
      expr = update_G_Rcpp(
        adj = adj,
        
        edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
        df = df_0 + n,
        df_0 = df_0,
        rate = diag(p) + U,
        n_edge = p,  # No. of single edge updates attempted in this MCMC update
        seed = sample.int(n = .Machine$integer.max, size = 1L),
        loc_bal = FALSE
      ), error = function (e) {
        warning(paste("`update_G_Rcpp` failed:", e))
        return (adj)
      }
    ))
  }
  
  
  p <- p_vec[p_ind]
  G_true <- G_true_list[[p_ind]][prop_ind, r, , ]
  data <- data_list[[p_ind]][prop_ind, r, , ]
  
  get_post_edge_prob <- function(n) {
    data_n <- data[seq_len(n), , drop = FALSE]
    U <- t(data_n) %*% data_n
    
    MCMC_step <- function(adj) update_G(
      adj = adj, edge_prob = 0.5,
      df_0 = 3.0, U = U, n = n
    )
    
    adj <- G_true
    for (s in 1:(n_iter %/% 10L)) adj <- MCMC_step(adj)  # Burn-in iterations
    adj_cum <- matrix(0L, nrow = p, ncol = p)
    
    for (s in 1:n_iter) {
      adj <- MCMC_step(adj)
      adj_cum <- adj_cum + adj
    }
    
    return (adj_cum / n_iter)
  }
  
  return (get_post_edge_prob(n_vec[n_ind]))
}


n_par <- n_rep * n_prop * n_n


run_task_par <- function(ind, p_ind) {
  n_ind <- (ind - 1L) %% n_n + 1L
  prop_ind <- ((ind - 1L) %/% n_n) %% n_prop + 1L
  r <- (ind - 1L) %/% (n_prop * n_n) + 1L
  return (run_task(prop_ind, r, n_ind, p_ind))
}


cl <- parallel::makeCluster(n_cores)
doSNOW::registerDoSNOW(cl)

for (p_ind in 1:n_p) {
  p <- p_vec[p_ind]
  print(paste("Working on p =", p))
  set.seed(1L)
  pb <- txtProgressBar(max = n_par, style = 3L)
  progress <- function (ind) setTxtProgressBar(pb, ind)
  
  result_list <- doRNG::`%dorng%`(obj = foreach::foreach(
    ind = 1:n_par, .options.snow = list(progress = progress)
  ), ex = run_task_par(ind, p_ind))
  
  close(pb)
  
  save(result_list, file = paste0("FPR_simul_p", p, "_result.RData"))
}

parallel::stopCluster(cl)



# Gather data in correct format.
inc_prob_list <- vector("list", n_p)
n_par <- n_rep * n_prop * n_n

for (p_ind in 1:n_p) {
  p <- p_vec[p_ind]
  load(paste0("FPR_simul_p", p, "_result.RData"))
  
  inc_prob_list[[p_ind]] <- array(
    data = NA_real_, dim = c(n_prop, n_rep, n_n, p, p)
  )
  
  for (ind in 1:n_par) {
    n_ind <- (ind - 1L) %% n_n + 1L
    prop_ind <- ((ind - 1L) %/% n_n) %% n_prop + 1L
    r <- (ind - 1L) %/% (n_prop * n_n) + 1L
    inc_prob_list[[p_ind]][prop_ind, r, n_ind, , ] <- result_list[[ind]]
  }
}



# Compute results for p = 100 using parallel computing.
print("Working on p = 100")
p <- 100L


run_task <- function (prop_ind, r, n_ind) {
  withr::with_makevars(
    new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
    code = Rcpp::sourceCpp(file = "src/ggm.cpp")
  )
  
  
  ## MCMC update for a graph
  
  # We now use the function `update_G_Rcpp` to define an R function that
  # provides an MCMC update for a graph in a Gaussian graphical model.
  
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
  update_G <- function(adj, edge_prob, df_0, U, n) {
    p <- nrow(adj)
    
    return (tryCatch(
      expr = update_G_Rcpp(
        adj = adj,
        
        edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
        df = df_0 + n,
        df_0 = df_0,
        rate = diag(p) + U,
        n_edge = p,  # No. of single edge updates attempted in this MCMC update
        seed = sample.int(n = .Machine$integer.max, size = 1L),
        loc_bal = FALSE
      ), error = function (e) {
        warning(paste("`update_G_Rcpp` failed:", e))
        return (adj)
      }
    ))
  }
  
  
  G_true <- G_true_arr[prop_ind, r, , ]
  data <- data_arr[prop_ind, r, , ]
  
  get_post_edge_prob <- function(n) {
    data_n <- data[seq_len(n), , drop = FALSE]
    U <- t(data_n) %*% data_n
    
    MCMC_step <- function(adj) update_G(
      adj = adj, edge_prob = 0.5,
      df_0 = 3.0, U = U, n = n
    )
    
    adj <- G_true
    for (s in 1:(n_iter %/% 10L)) adj <- MCMC_step(adj)  # Burn-in iterations
    adj_cum <- matrix(0L, nrow = p, ncol = p)
    
    for (s in 1:n_iter) {
      adj <- MCMC_step(adj)
      adj_cum <- adj_cum + adj
    }
    
    return (adj_cum / n_iter)
  }
  
  return (get_post_edge_prob(n_vec[n_ind]))
}


n_par <- n_rep * n_prop * n_n


run_task_par <- function(ind) {
  n_ind <- (ind - 1L) %% n_n + 1L
  prop_ind <- ((ind - 1L) %/% n_n) %% n_prop + 1L
  r <- (ind - 1L) %/% (n_prop * n_n) + 1L
  return (run_task(prop_ind, r, n_ind))
}


cl <- parallel::makeCluster(n_cores)
doSNOW::registerDoSNOW(cl)
set.seed(1L)
pb <- txtProgressBar(max = n_par, style = 3L)
progress <- function (ind) setTxtProgressBar(pb, ind)

result_list <- doRNG::`%dorng%`(obj = foreach::foreach(
  ind = 1:n_par, .options.snow = list(progress = progress)
), ex = run_task_par(ind))

close(pb)
parallel::stopCluster(cl)



# Gather data in correct format.
p_vec <- c(20L, 40L, p)
n_p <- length(p_vec)
G_true_list[[n_p]] <- G_true_arr

inc_prob_list[[n_p]] <- array(
  data = NA_real_, dim = c(n_prop, n_rep, n_n, p, p)
)


for (ind in 1:n_par) {
  n_ind <- (ind - 1L) %% n_n + 1L
  prop_ind <- ((ind - 1L) %/% n_n) %% n_prop + 1L
  r <- (ind - 1L) %/% (n_prop * n_n) + 1L
  inc_prob_list[[n_p]][prop_ind, r, n_ind, , ] <- result_list[[ind]]
}



# Plot the results.
cols <- rev(viridisLite::viridis(n_prop, alpha = 0.9))
lwds <- rep(2L, n_prop)
ltys <- c(2L, 1L, 2L, 1L)
pdf(file = "simul_ROC.pdf", width = 6.5, height = 8)
par(mfcol = c(n_n - 1L, n_p), mai = rep(0, 4L), omi = rep(1, 4L))

for (p_ind in 1:n_p) for (n_ind in 2:n_n) for (prop_ind in 2:n_prop) {
  p <- p_vec[p_ind]
  ind <- upper.tri(diag(p))
  inc_prob_l <- list(numeric(0L), numeric(0))
  
  for (r in 1:n_rep) for (tmp_ind in 1:2) inc_prob_l[[tmp_ind]] <- c(
    inc_prob_l[[tmp_ind]],
    inc_prob_list[[p_ind]][prop_ind, r, n_ind, , ][
      ind & G_true_list[[p_ind]][prop_ind, r, , ] == tmp_ind - 1L
    ]
  )
  
  inc_prob <- c(
    1, sort(unique(inc_prob_l[[1L]], inc_prob_l[[2L]]), decreasing = TRUE), 0
  )
  
  n_inc_prob <- length(inc_prob)
  x <- numeric(n_inc_prob)
  y <- numeric(n_inc_prob)
  
  for (i in 1:n_inc_prob) {
    x[i] <- mean(inc_prob_l[[1L]] > inc_prob[i])
    y[i] <- mean(inc_prob_l[[2L]] > inc_prob[i])
  }
  
  if (prop_ind == 2L) {
    plot(
      x = 0, type = "n", ylab = "True positive rate", xaxs = "i", yaxs = "i",
      yaxt = if (p_ind == 1L) "s" else "n",
      xaxt = if (n_ind == n_n) "s" else "n", xaxp = c(0.2, 0.8, 3L),
      yaxp = c(0.2, 0.8, 3L), xlim = 0:1, ylim = 0:1, asp = 1
    )
    
    abline(a = 0, b = 1, col = "darkgrey")
  }
  
  if (n_ind == 2L) mtext(
    text = paste("p =", p), side = 3L, line = 1L, outer = TRUE,
    at = 2 * p_ind / (2 * n_p + 1) * 1.16 - 0.15
  )
  
  if (n_ind == n_n) mtext(
    text = "False positive rate", side = 1L, line = 3L,
    outer = TRUE, at = 2 * p_ind / (2 * n_p + 1) * 1.16 - 0.17
  )
  
  if (p_ind == n_p) mtext(
    text = paste("n =", n_vec[n_ind], "\nTrue positive rate"), side = 2L,
    line = 4L, outer = TRUE,
    at = 2 * (n_n + 1L - n_ind) / (2 * n_n + 1) * 1.37 - 0.12
  )
  
  
  ROC <- stepfun(x = x, y = c(y, 1))
  
  curve(
    expr = ROC, from = 0, to = 1, lwd = lwds[prop_ind], lty = ltys[prop_ind],
    col = cols[prop_ind], add = TRUE
  )
  
  points(
    x = mean(inc_prob_l[[1L]] > 0.5), y = mean(inc_prob_l[[2L]] > 0.5),
    pch = 16L, col = cols[prop_ind], cex = 1.5
  )
  
  if (p_ind == n_p & prop_ind == 2L & n_ind == n_n) legend(
    x = "bottomright",
    legend = c("No edges", paste0(100 * prop_vec[-1], "% of edges"))[-1],
    lty = ltys[-1], lwd = lwds[-1], col = cols[-1]
  )
}

dev.off()

