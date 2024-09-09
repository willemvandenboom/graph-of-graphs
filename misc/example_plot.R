### Graph of graphs: Example

# Code for the example on how the graph of graphs differs from standard Gaussian
# graphical models (GGMs) in the paper "Graph of Graphs: From Nodes to
# Supernodes in Graphical Models" by Maria De Iorio, Willem van den Boom,
# Alexandros Beskos, Ajay Jasra and Andrea Cremaschi.


source("graph_of_graphs_MCMC.R")
n <- 1e3L
supernode_sizes <- c(3L, 5L, 5L, 6L)
n_iter <- 2e2L

superedge <- TRUE
superedge_strength <- -.2

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
  prec_mat[i1, i2] <- 0.1 / sqrt(supernode_sizes[k] - 1L)
  
  prec_mat[i2, i1] <- prec_mat[i1, i2]
}

if (superedge) {
  i1 <- center_ind[K - 1L]
  i2 <- center_ind[K]
  prec_mat[i1, i2] <- superedge_strength
  prec_mat[i2, i1] <- prec_mat[i1, i2]
}

set.seed(1L)

# Make the first supernode complete
prec_mat[2, 3] <- prec_mat[1, 2]
prec_mat[3, 2] <- prec_mat[2, 3]

# Add more "superedges"
prec_mat[i2 + 1:3, i1] <- prec_mat[i1, i2]
prec_mat[i1, i2 + 1:3] <- prec_mat[i1, i2]

prec_mat[1, p] <- superedge_strength
prec_mat[p, 1] <- prec_mat[1, p]


for (i in 1:3) {
  prec_mat[supernode_sizes[1] + i, sum(supernode_sizes[1:2]) + i] <- -.6
  prec_mat[sum(supernode_sizes[1:2]) + i, + supernode_sizes[1] + i] <- -.6
}

X <- mvtnorm::rmvnorm(n = n, sigma = cov2cor(solve(prec_mat)))


# Set tessellation to the truth
assignment_nested <- rep(x = 1:K, times = supernode_sizes)



print("Plotting the graph sphere...")


# Estimate the supergraph with the tessellation fixed.


# Function that runs MCMC and generates output for the supergraph only
run_MCMC_G_only <- function (
    X, assignment, n_iter, burnin = NULL
) {
  if (is.null(burnin)) burnin <- n_iter %/% 10L
  n <- nrow(X)
  n_centers <- max(assignment)
  G_MCMC <- array(data = NA_integer_, dim = c(n_iter, n_centers, n_centers))
  Omega_MCMC <- array(data = NA_real_, dim = c(n_iter, n_centers, n_centers))
  super_edge_mat <- matrix(data = 0L, nrow = n_centers, ncol = n_centers)
  U <- compute_U_unaugmented(X = X, assignment = assignment)
  
  # Supergraph: initialize at empty
  G <- matrix(0L, nrow = n_centers, ncol = n_centers)
  
  pb <- txtProgressBar(max = burnin + n_iter, style = 3)
  
  for (s in 1:(burnin + n_iter)) {
    G <- update_G_WWA(G, 0.5, 3, U, n, 1L)
    
    if (s > burnin) {
      G_MCMC[s - burnin, , ] <- G
      Omega_MCMC[s - burnin, , ] <- rgwish_Rcpp(G, 3 + n, diag(n_centers) + U, s)
      super_edge_mat <- super_edge_mat + G
    }
    
    setTxtProgressBar(pb, s)
  }
  
  close(pb)
  
  return(list(
    G_MCMC = G_MCMC, Omega_MCMC = Omega_MCMC,
    superedge_mat = super_edge_mat / n_iter
  ))
}



print("Running MCMC on the supergraph with the tessellation fixed...")
set.seed(1L)
res_G <- run_MCMC_G_only(X = X, assignment = assignment_nested, n_iter = 1e5L)
supergraph_adj <- res_G$superedge_mat > 0.5

print("Point estimate of the supergraph with 0.5 edge inclusion threshold:")

print(igraph::graph_from_adjacency_matrix(
  adjmatrix = supergraph_adj, mode = "undirected", diag = FALSE
))

res_G$superedge_mat



mod_sizes <- supernode_sizes


p_k_sub <- sort(table(assignment_nested))
p_sub <- sum(p_k_sub)
reorder_supernodes <- as.integer(names(p_k_sub))
K_nested <- 4L

super_layout <- matrix(data = 0, nrow = K_nested, ncol = 2L)
super_layout[2:3, 1] <- 1
super_layout[3:4, 2] <- 1


G_list <- vector(mode = "list", length = K_nested)
layout <- matrix(NA_real_, nrow = p_sub, ncol = 2L)

# Compute MAP (maximum a posteriori) estimates of the trees using `log_w_matrix`
# These estimates do not depend on the coarsening parameter `coars_pow`.
log_w_matrix <- compute_log_w_matrix(X = X, coars_pow = 1)

for (k in 1:K_nested) {
  ind <- which(assignment_nested == reorder_supernodes[k])
  
  G_list[[k]] <- igraph::mst(igraph::graph_from_adjacency_matrix(
    adjmatrix = -log_w_matrix[ind, ind, drop = FALSE], mode = "undirected",
    weighted = TRUE, diag = FALSE
  ))
  
  igraph::V(G_list[[k]])$name <- as.character(ind)
  
  # Permute the nodes using depth-first search for pretty plotting.
  G_list[[k]] <- igraph::permute(graph = G_list[[k]], permutation = order(
    igraph::dfs(
      graph = G_list[[k]], root = which.min(igraph::degree(G_list[[k]]))
    )$order
  ))
  
  layout[sum(p_k_sub[seq_len(k - 1L)]) + 1:p_k_sub[k], ] <- t(
    5 * super_layout[k, ] + t(igraph::layout_as_tree(
      graph = G_list[[k]], circular = TRUE
    ))
  )[order(as.integer(names(igraph::V(G_list[[k]])))), ]
  
  if (k == 4L) layout[sum(p_k_sub[seq_len(k - 1L)]) + 1:p_k_sub[k], ] <- t(
    5 * super_layout[k, ] + t(igraph::layout_as_star(
      graph = G_list[[k]]
    ))
  )[order(as.integer(names(igraph::V(G_list[[k]])))), ]
}

G_block <- igraph:::union.igraph(G_list, byname = TRUE)

cols <- RColorBrewer::brewer.pal(
  n = 6L, name = "Set1"
)[c(3L, 1L, 2L, 4L, 5L, 6L)]

mod_ind <- rep(1:length(mod_sizes), mod_sizes)

# Ensure that the vertices are labeled consecutively.
G_block <- igraph::permute(
  G_block, order(as.integer(names(igraph::V(G_block))))
)

tmp <- list(membership = assignment_nested[as.integer(igraph::V(G_block)$name)])
class(tmp) <- "communities"



library(igraph)


# Modified version of `igraph::plot.igraph` to visualize a graph sphere
plot.igraph2 <- function (x, supergraph_adj, axes = FALSE, add = FALSE, xlim = c(-1, 1), ylim = c(-1, 
                                                                                                  1), mark.groups = list(), mark.shape = 1/2, mark.col = rainbow(length(mark.groups), 
                                                                                                                                                                 alpha = 0.3), mark.border = rainbow(length(mark.groups), 
                                                                                                                                                                                                     alpha = 1), mark.expand = 15, loop.size = 1, ...) 
{
  graph <- x
  igraph:::ensure_igraph(graph)
  vc <- vcount(graph)
  params <- igraph:::i.parse.plot.params(graph, list(...))
  vertex.size <- 1/200 * params("vertex", "size")
  label.family <- params("vertex", "label.family")
  label.font <- params("vertex", "label.font")
  label.cex <- params("vertex", "label.cex")
  label.degree <- params("vertex", "label.degree")
  label.color <- params("vertex", "label.color")
  label.dist <- params("vertex", "label.dist")
  labels <- params("vertex", "label")
  shape <- igraph:::igraph.check.shapes(params("vertex", "shape"))
  edge.color <- params("edge", "color")
  edge.width <- params("edge", "width")
  edge.lty <- params("edge", "lty")
  arrow.mode <- params("edge", "arrow.mode")
  edge.labels <- params("edge", "label")
  loop.angle <- params("edge", "loop.angle")
  edge.label.font <- params("edge", "label.font")
  edge.label.family <- params("edge", "label.family")
  edge.label.cex <- params("edge", "label.cex")
  edge.label.color <- params("edge", "label.color")
  elab.x <- params("edge", "label.x")
  elab.y <- params("edge", "label.y")
  arrow.size <- params("edge", "arrow.size")[1]
  arrow.width <- params("edge", "arrow.width")[1]
  curved <- params("edge", "curved")
  if (is.function(curved)) {
    curved <- curved(graph)
  }
  layout <- igraph:::i.postprocess.layout(params("plot", "layout"))
  margin <- params("plot", "margin")
  margin <- rep(margin, length.out = 4)
  rescale <- params("plot", "rescale")
  asp <- params("plot", "asp")
  frame.plot <- params("plot", "frame.plot")
  main <- params("plot", "main")
  sub <- params("plot", "sub")
  xlab <- params("plot", "xlab")
  ylab <- params("plot", "ylab")
  palette <- params("plot", "palette")
  if (!is.null(palette)) {
    old_palette <- palette(palette)
    on.exit(palette(old_palette), add = TRUE)
  }
  arrow.mode <- igraph:::i.get.arrow.mode(graph, arrow.mode)
  maxv <- max(vertex.size)
  if (vc > 0 && rescale) {
    layout <- norm_coords(layout, -1, 1, -1, 1)
    xlim <- c(xlim[1] - margin[2] - maxv, xlim[2] + margin[4] + 
                maxv)
    ylim <- c(ylim[1] - margin[1] - maxv, ylim[2] + margin[3] + 
                maxv)
  }
  if (!add) {
    plot(0, 0, type = "n", xlab = xlab, ylab = ylab, xlim = xlim, 
         ylim = ylim, axes = axes, frame.plot = ifelse(is.null(frame.plot), 
                                                       axes, frame.plot), asp = asp, main = main, sub = sub)
  }
  if (!is.list(mark.groups) && is.numeric(mark.groups)) {
    mark.groups <- list(mark.groups)
  }
  if (inherits(mark.groups, "communities")) {
    mark.groups <- communities(mark.groups)
  }
  mark.shape <- rep(mark.shape, length.out = length(mark.groups))
  mark.border <- rep(mark.border, length.out = length(mark.groups))
  mark.col <- rep(mark.col, length.out = length(mark.groups))
  mark.expand <- rep(mark.expand, length.out = length(mark.groups))
  K <- length(mark.groups)
  supernode_layout <- matrix(data = NA_real_, nrow = K, ncol = 2L)
  
  for (g in seq_along(mark.groups)) {
    .members <- mark.groups[[g]]
    v <- V(graph)[.members]
    supernode_layout[g, ] <- colMeans(layout[v, , drop = FALSE])
  }
  
  s_ind <- as.integer(names(communities(tmp)))
  lwd <- 2
  
  if (K > 1L) for (k1 in 1:(K - 1L)) for (k2 in (k1 + 1L):K) if (
    supergraph_adj[s_ind[k1], s_ind[k2]] | supergraph_adj[s_ind[k2], s_ind[k1]]
  ) segments(
    x0 = supernode_layout[k1, 1], y0 = supernode_layout[k1, 2],
    x1 = supernode_layout[k2, 1], y1 = supernode_layout[k2, 2], lwd = lwd,
    col = mark.border[1]
  )
  
  points(x = supernode_layout[, 1], y = supernode_layout[, 2], col = mark.border[1], pch = 15L, cex = 2)
  
  igraph.polygon2 <- function (points, vertex.size = 15/200, expand.by = 15/200, shape = 1/2, 
                               col = "#ff000033", border = NA) {
    by <- expand.by
    pp <- rbind(points, cbind(points[, 1] - vertex.size - by, 
                              points[, 2]), cbind(points[, 1] + vertex.size + by, points[, 
                                                                                         2]), cbind(points[, 1], points[, 2] - vertex.size - by), 
                cbind(points[, 1], points[, 2] + vertex.size + by))
    cl <- convex_hull(pp)
    xspline(cl$rescoords, shape = shape, open = FALSE, col = col, 
            border = border, lwd = lwd)
  }
  
  for (g in seq_along(mark.groups)) {
    .members <- mark.groups[[g]]
    v <- V(graph)[.members]
    if (length(vertex.size) == 1) {
      vs <- vertex.size
    }
    else {
      vs <- rep(vertex.size, length.out = vcount(graph))[v]
    }
    igraph.polygon2(layout[v, , drop = FALSE], vertex.size = vs, 
                    expand.by = mark.expand[g]/200, shape = mark.shape[g], 
                    col = mark.col[g], border = mark.border[g])
  }
  
  el <- as_edgelist(graph, names = FALSE)
  loops.e <- which(el[, 1] == el[, 2])
  nonloops.e <- which(el[, 1] != el[, 2])
  loops.v <- el[, 1][loops.e]
  loop.labels <- edge.labels[loops.e]
  loop.labx <- if (is.null(elab.x)) {
    rep(NA, length(loops.e))
  }
  else {
    elab.x[loops.e]
  }
  loop.laby <- if (is.null(elab.y)) {
    rep(NA, length(loops.e))
  }
  else {
    elab.y[loops.e]
  }
  edge.labels <- edge.labels[nonloops.e]
  elab.x <- if (is.null(elab.x)) 
    NULL
  else elab.x[nonloops.e]
  elab.y <- if (is.null(elab.y)) 
    NULL
  else elab.y[nonloops.e]
  el <- el[nonloops.e, , drop = FALSE]
  edge.coords <- matrix(0, nrow = nrow(el), ncol = 4)
  edge.coords[, 1] <- layout[, 1][el[, 1]]
  edge.coords[, 2] <- layout[, 2][el[, 1]]
  edge.coords[, 3] <- layout[, 1][el[, 2]]
  edge.coords[, 4] <- layout[, 2][el[, 2]]
  if (length(unique(shape)) == 1) {
    ec <- igraph:::.igraph.shapes[[shape[1]]]$clip(edge.coords, el, 
                                                   params = params, end = "both")
  }
  else {
    shape <- rep(shape, length.out = vcount(graph))
    ec <- edge.coords
    ec[, 1:2] <- t(sapply(seq(length.out = nrow(el)), function(x) {
      igraph:::.igraph.shapes[[shape[el[x, 1]]]]$clip(edge.coords[x, 
                                                                  , drop = FALSE], el[x, , drop = FALSE], params = params, 
                                                      end = "from")
    }))
    ec[, 3:4] <- t(sapply(seq(length.out = nrow(el)), function(x) {
      igraph:::.igraph.shapes[[shape[el[x, 2]]]]$clip(edge.coords[x, 
                                                                  , drop = FALSE], el[x, , drop = FALSE], params = params, 
                                                      end = "to")
    }))
  }
  x0 <- ec[, 1]
  y0 <- ec[, 2]
  x1 <- ec[, 3]
  y1 <- ec[, 4]
  if (length(loops.e) > 0) {
    ec <- edge.color
    if (length(ec) > 1) {
      ec <- ec[loops.e]
    }
    point.on.cubic.bezier <- function(cp, t) {
      c <- 3 * (cp[2, ] - cp[1, ])
      b <- 3 * (cp[3, ] - cp[2, ]) - c
      a <- cp[4, ] - cp[1, ] - c - b
      t2 <- t * t
      t3 <- t * t * t
      a * t3 + b * t2 + c * t + cp[1, ]
    }
    compute.bezier <- function(cp, points) {
      dt <- seq(0, 1, by = 1/(points - 1))
      sapply(dt, function(t) point.on.cubic.bezier(cp, 
                                                   t))
    }
    plot.bezier <- function(cp, points, color, width, arr, 
                            lty, arrow.size, arr.w) {
      p <- compute.bezier(cp, points)
      polygon(p[1, ], p[2, ], border = color, lwd = width, 
              lty = lty)
      if (arr == 1 || arr == 3) {
        igraph:::igraph.Arrows(p[1, ncol(p) - 1], p[2, ncol(p) - 
                                                      1], p[1, ncol(p)], p[2, ncol(p)], sh.col = color, 
                               h.col = color, size = arrow.size, sh.lwd = width, 
                               h.lwd = width, open = FALSE, code = 2, width = arr.w)
      }
      if (arr == 2 || arr == 3) {
        igraph:::igraph.Arrows(p[1, 2], p[2, 2], p[1, 1], p[2, 
                                                            1], sh.col = color, h.col = color, size = arrow.size, 
                               sh.lwd = width, h.lwd = width, open = FALSE, 
                               code = 2, width = arr.w)
      }
    }
    loop <- function(x0, y0, cx = x0, cy = y0, color, angle = 0, 
                     label = NA, width = 1, arr = 2, lty = 1, arrow.size = arrow.size, 
                     arr.w = arr.w, lab.x, lab.y, loopSize = loop.size) {
      rad <- angle
      center <- c(cx, cy)
      cp <- matrix(c(x0, y0, x0 + 0.4 * loopSize, y0 + 
                       0.2 * loopSize, x0 + 0.4 * loopSize, y0 - 0.2 * 
                       loopSize, x0, y0), ncol = 2, byrow = TRUE)
      phi <- atan2(cp[, 2] - center[2], cp[, 1] - center[1])
      r <- sqrt((cp[, 1] - center[1])^2 + (cp[, 2] - center[2])^2)
      phi <- phi + rad
      cp[, 1] <- cx + r * cos(phi)
      cp[, 2] <- cy + r * sin(phi)
      if (is.na(width)) {
        width <- 1
      }
      plot.bezier(cp, 50, color, width, arr = arr, lty = lty, 
                  arrow.size = arrow.size, arr.w = arr.w)
      if (is.language(label) || !is.na(label)) {
        lx <- x0 + 0.3
        ly <- y0
        phi <- atan2(ly - center[2], lx - center[1])
        r <- sqrt((lx - center[1])^2 + (ly - center[2])^2)
        phi <- phi + rad
        lx <- cx + r * cos(phi)
        ly <- cy + r * sin(phi)
        if (!is.na(lab.x)) {
          lx <- lab.x
        }
        if (!is.na(lab.y)) {
          ly <- lab.y
        }
        text(lx, ly, label, col = edge.label.color, font = edge.label.font, 
             family = edge.label.family, cex = edge.label.cex)
      }
    }
    ec <- edge.color
    if (length(ec) > 1) {
      ec <- ec[loops.e]
    }
    vs <- vertex.size
    if (length(vertex.size) > 1) {
      vs <- vs[loops.v]
    }
    ew <- edge.width
    if (length(edge.width) > 1) {
      ew <- ew[loops.e]
    }
    la <- loop.angle
    if (length(loop.angle) > 1) {
      la <- la[loops.e]
    }
    lty <- edge.lty
    if (length(edge.lty) > 1) {
      lty <- lty[loops.e]
    }
    arr <- arrow.mode
    if (length(arrow.mode) > 1) {
      arr <- arrow.mode[loops.e]
    }
    asize <- arrow.size
    if (length(arrow.size) > 1) {
      asize <- arrow.size[loops.e]
    }
    xx0 <- layout[loops.v, 1] + cos(la) * vs
    yy0 <- layout[loops.v, 2] - sin(la) * vs
    mapply(loop, xx0, yy0, color = ec, angle = -la, label = loop.labels, 
           lty = lty, width = ew, arr = arr, arrow.size = asize, 
           arr.w = arrow.width, lab.x = loop.labx, lab.y = loop.laby)
  }
  if (length(x0) != 0) {
    if (length(edge.color) > 1) {
      edge.color <- edge.color[nonloops.e]
    }
    if (length(edge.width) > 1) {
      edge.width <- edge.width[nonloops.e]
    }
    if (length(edge.lty) > 1) {
      edge.lty <- edge.lty[nonloops.e]
    }
    if (length(arrow.mode) > 1) {
      arrow.mode <- arrow.mode[nonloops.e]
    }
    if (length(arrow.size) > 1) {
      arrow.size <- arrow.size[nonloops.e]
    }
    if (length(curved) > 1) {
      curved <- curved[nonloops.e]
    }
    if (length(unique(arrow.mode)) == 1) {
      lc <- igraph:::igraph.Arrows(x0, y0, x1, y1, h.col = edge.color, 
                                   sh.col = edge.color, sh.lwd = edge.width, h.lwd = 1, 
                                   open = FALSE, code = arrow.mode[1], sh.lty = edge.lty, 
                                   h.lty = 1, size = arrow.size, width = arrow.width, 
                                   curved = curved)
      lc.x <- lc$lab.x
      lc.y <- lc$lab.y
    }
    else {
      curved <- rep(curved, length.out = ecount(graph))[nonloops.e]
      lc.x <- lc.y <- numeric(length(curved))
      for (code in 0:3) {
        valid <- arrow.mode == code
        if (!any(valid)) {
          next
        }
        ec <- edge.color
        if (length(ec) > 1) {
          ec <- ec[valid]
        }
        ew <- edge.width
        if (length(ew) > 1) {
          ew <- ew[valid]
        }
        el <- edge.lty
        if (length(el) > 1) {
          el <- el[valid]
        }
        lc <- igraph:::igraph.Arrows(x0[valid], y0[valid], x1[valid], 
                                     y1[valid], code = code, sh.col = ec, h.col = ec, 
                                     sh.lwd = ew, h.lwd = 1, h.lty = 1, sh.lty = el, 
                                     open = FALSE, size = arrow.size, width = arrow.width, 
                                     curved = curved[valid])
        lc.x[valid] <- lc$lab.x
        lc.y[valid] <- lc$lab.y
      }
    }
    if (!is.null(elab.x)) {
      lc.x <- ifelse(is.na(elab.x), lc.x, elab.x)
    }
    if (!is.null(elab.y)) {
      lc.y <- ifelse(is.na(elab.y), lc.y, elab.y)
    }
    text(lc.x, lc.y, labels = edge.labels, col = edge.label.color, 
         family = edge.label.family, font = edge.label.font, 
         cex = edge.label.cex)
  }
  rm(x0, y0, x1, y1)
  if (vc > 0) {
    if (length(unique(shape)) == 1) {
      igraph:::.igraph.shapes[[shape[1]]]$plot(layout, params = params)
    }
    else {
      sapply(seq(length.out = vcount(graph)), function(x) {
        igraph:::.igraph.shapes[[shape[x]]]$plot(layout[x, , drop = FALSE], 
                                                 v = x, params = params)
      })
    }
  }
  old_xpd <- par(xpd = TRUE)
  on.exit(par(old_xpd), add = TRUE)
  x <- layout[, 1] + label.dist * cos(-label.degree) * (vertex.size + 
                                                          6 * 8 * log10(2))/200
  y <- layout[, 2] + label.dist * sin(-label.degree) * (vertex.size + 
                                                          6 * 8 * log10(2))/200
  if (vc > 0) {
    if (length(label.family) == 1) {
      text(x, y, labels = labels, col = label.color, family = label.family, 
           font = label.font, cex = label.cex)
    }
    else {
      if1 <- function(vect, idx) if (length(vect) == 1) 
        vect
      else vect[idx]
      sapply(seq_len(vcount(graph)), function(v) {
        text(x[v], y[v], labels = if1(labels, v), col = if1(label.color, 
                                                            v), family = if1(label.family, v), font = if1(label.font, 
                                                                                                          v), cex = if1(label.cex, v))
      })
    }
  }
  rm(x, y)
  invisible(NULL)
}




pdf_wh <- 4
vertex.size <- 5
pdf("GGM_to_graph_of_graphs.pdf", width = 2 * pdf_wh, height = pdf_wh)
par(mar = c(0, 0, 3, 0), mfrow = 1:2)

plot.igraph2(
  x = igraph::graph_from_adjacency_matrix(
    adjmatrix = prec_mat != 0, mode = "undirected", diag = FALSE
  ),
  supergraph_adj = matrix(data = FALSE, nrow = K_nested, ncol = K_nested),
  mark.groups = tmp, mark.col = "white", mark.border = "white",
  mark.expand = 30, vertex.size = vertex.size,
  vertex.color = cols[mod_ind[as.integer(V(G_block)$name)]],
  vertex.label = "", edge.color = "black", layout = layout[, 2:1],
  xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), main = "Latent graph"
)

plot.igraph2(
  x = G_block, supergraph_adj = supergraph_adj, mark.groups = tmp,
  mark.col = scales::alpha("white", 0.5), mark.border = "darkgray",
  mark.expand = 30, vertex.size = vertex.size,
  vertex.color = cols[mod_ind[as.integer(V(G_block)$name)]],
  vertex.label = "", edge.color = "black", layout = layout[, 2:1],
  xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), main = "Graph of graphs"
)

dev.off()
