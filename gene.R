### Graph sphere: Application to gene expression data

# Code for the application with gene expression data in the paper "Graph Sphere:
# From Nodes to Supernodes in Graphical Models" by Willem van den Boom, Maria De
# Iorio, Alexandros Beskos and Ajay Jasra.



print("Installing required packages...")
pkgs <- c("BiocManager", "ggplot2", "glasso", "httr", "igraph", "remotes")
pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(pkgs) > 0L) install.packages(
  pkgs = pkgs, dependencies = TRUE, repos = "https://cloud.r-project.org"
)

pkgs <- c(
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "readxl", "RTCGA",
  "RTCGA.mRNA", "STRINGdb"
)

pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(pkgs) > 0L) {
  BiocManager::install(pkgs = pkgs, dependencies = TRUE, update = FALSE)
}



print("Downloading the gene expression data...")


# The following function to read and transpose a CSV file is taken from
# https://gist.github.com/jonocarroll/b17ce021b0637a31f584ed08a1fbe733#file-read-tscv-r-L1-L20.
# Based on https://stackoverflow.com/a/17289991/4168169
read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}


# Additional file 2 of Zhang (2018, doi:10.1186/s12918-018-0530-9)
mods <- read.tcsv(
  file = "https://static-content.springer.com/esm/art%3A10.1186%2Fs12918-018-0530-9/MediaObjects/12918_2018_530_MOESM2_ESM.csv",
  sep = ";"
)

# Select genes based on Table 2 in Zhang (2018, doi:10.1186/s12918-018-0530-9):
mods_selected <- gsub(" ", "", mods$modulenumberinfilename[
  mods$Modulenumberinpaper %in% c(5L, 6L, 14L, 27L, 36L, 39L)
])

# Additional file 1 of Zhang (2018, doi:10.1186/s12918-018-0530-9)
httr::GET(
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs12918-018-0530-9/MediaObjects/12918_2018_530_MOESM1_ESM.xlsx",
  httr::write_disk(temp_file <- tempfile(fileext = ".xlsx"))
)

n_mods <- length(mods_selected)
genes_selected <- list()
all_genes <- character(0)

for (i in 1:n_mods) {
  genes_selected[[i]] <- readxl::read_xlsx(
    path = temp_file, sheet = mods_selected[i]
  )[[2]]
  
  all_genes <- c(all_genes, genes_selected[[i]])
}

unlink(temp_file)

mod_sizes <- integer(n_mods)  # Sizes of the modules
for (i in 1:n_mods) mod_sizes[i] <- length(genes_selected[[i]])

tmp <- RTCGA::expressionsTCGA(
  RTCGA.mRNA::BRCA.mRNA, RTCGA.mRNA::OV.mRNA, extract.cols = all_genes
)

tmp$dataset <- gsub(pattern = ".mRNA", replacement = "",  tmp$dataset)
tmp$dataset <- gsub(pattern = "RTCGA.mRNA::", replacement = "",  tmp$dataset)
tmp$dataset <- gsub(pattern = "RTCGA::", replacement = "",  tmp$dataset)

X <- as.matrix(subset(
  tmp[complete.cases(tmp), c("dataset", all_genes)][, -1], tmp$dataset == "OV"
))


# Quantile-normalize the data to marginally follow a standard Gaussian
# distribution.
quantile_normalize <- function (x) {
  qnorm(rank(x = x, ties.method = "average") / (length(x) + 1L))
}


X <- apply(X = X, MARGIN = 2L, FUN = quantile_normalize)



source("graph_sphere_MCMC.R")
p <- ncol(X)
log_cohesion = dnbinom(x = 0:(p - 1L), size = 2, prob = 1 / 6, log = TRUE)
n_iter <- 1e5L

set.seed(1L)
print("Running MCMC with coarsening of the likelihood...")

res_coars <- run_MCMC(
  X = X, nested = FALSE, n_iter = n_iter, log_cohesion = log_cohesion
)

print("Running nested MCMC...")

res_nested <- run_MCMC(
  X = X, nested = TRUE, n_iter = n_iter, log_cohesion = log_cohesion
)



print("Creating trace plots...")
pdf("trace_MCMC.pdf", width = 8L, height = 5L)
options(scipen=10000)
par(mfcol = c(2L, 2L), mar = c(4, 5, 3, 2), xaxs = "i")

for (i in 1:2) {
  res <- if (i == 1L) res_coars else res_nested
  res$K_MCMC <- integer(n_iter)
  for (s in 1:n_iter) res$K_MCMC[s] <- NCOL(res$G_MCMC[[s]])
  xlab <- if (i == 1L) "MCMC iteration" else "Outer MCMC iteration"
  
  plot(
    x = res$K_MCMC, type = "l", xlab = xlab,
    ylab = substitute(paste("Number of supernodes ", italic("K"))),
    main = if(i == 1L) "Coarsened likelihood" else "Nested MCMC"
  )
  
  plot(
    x = res$log_post_MCMC, type = "l", xlab = xlab,
    ylab = if (i == 1L) {
      expression("log{"*pi^(zeta)*(italic(T)*", "*italic(G)*", "*italic(Y))*"}")
    } else expression("log{"*italic(p)^(zeta)*(italic(T))*"}")
  )
}

dev.off()



print("Plotting the posterior coclustering  and superedge probabilities...")


# Function to reorder and plot similarity function
# The reordering uses hierarchical clustering.
get_order <- function (res) {
  sim_mat <- res$post_sim
  ind_reorder <- integer(p)
  
  # Reorder each module separately
  for (i in 1:n_mods) {
    start_ind <- sum(mod_sizes[seq_len(i - 1L)])
    ind_block <- start_ind + 1:mod_sizes[i]
    
    ind_reorder[ind_block] <- start_ind + heatmap(
      sim_mat[ind_block, ind_block]
    )$rowInd
  }
  
  return (ind_reorder)
}


plot_sim_mat <- function (
    res, ind_reorder, sim_mat = NULL,
    label = "Posterior coclustering probability"
) {
  # Reorder the modules from small to large.
  block_order <- order(mod_sizes)
  ind_order <- integer(p)
  
  for (i in 1:n_mods) {
    block_size <- mod_sizes[block_order[i]]
    
    ind_order[sum(mod_sizes[block_order][seq_len(i - 1L)]) + 1:block_size] <-
      ind_reorder[sum(mod_sizes[seq_len(block_order[i] - 1L)]) + 1:block_size]
  }
  
  if (is.null(sim_mat)) sim_mat <- res$post_sim
  mat <- sim_mat[ind_order, ind_order]
  p <- NCOL(mat)
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
    colorkey = list(space = "bottom", title = label),
    panel = function (...) {
      lattice::panel.levelplot(...)
      
      for (ind in cumsum(mod_sizes[block_order][-n_mods])) {
        lattice::panel.abline(h = p - ind + 0.5, col = "red", lty = 2L)
        lattice::panel.abline(v = ind + 0.5, col = "red", lty = 2L)
      }
    }
  )
}


ind_reorder <- get_order(res_nested)
library(latticeExtra)
pdf("gene_post_sim.pdf", width = 8L, height = 5L)

print(c(
  "Coarsened likelihood" = plot_sim_mat(
    res = res_coars, ind_reorder = ind_reorder
  ), "Nested MCMC" = plot_sim_mat(res = res_nested, ind_reorder = ind_reorder),
  layout = c(2L, 1L), merge.legends = FALSE
))

dev.off()
pdf("gene_superedge.pdf", width = 8L, height = 5L)

print(c(
  "Coarsened likelihood" = plot_sim_mat(
    res = res_coars, ind_reorder = ind_reorder,
    sim_mat = res_coars$superedge_mat,
    label = "Posterior superedge probability"
  ), "Nested MCMC" = plot_sim_mat(
    res = res_nested, ind_reorder = ind_reorder,
    sim_mat = res_nested$superedge_mat
  ),
  layout = c(2L, 1L), merge.legends = FALSE
))

dev.off()



print("Computing point estimate of the tessellation...")
remotes::install_github("sarawade/mcclust.ext")

X_cor <- abs(cor(X))
res_coars$assignment_MCMC <- matrix(data = NA_integer_, nrow = n_iter, ncol = p)

res_nested$assignment_MCMC <- matrix(
  data = NA_integer_, nrow = n_iter, ncol = p
)

for (s in 1:n_iter) {
  res_coars$assignment_MCMC[s, ] <- centers2tessellation(
    X_cor,  res_coars$centers_MCMC[s, ]
  )
  
  res_nested$assignment_MCMC[s, ] <- centers2tessellation(
    X_cor,  res_nested$centers_MCMC[s, ]
  )
}

assignment_coars <- mcclust.ext::minVI(
  psm = mcclust::comp.psm(res_coars$assignment),
  cls.draw = res_coars$assignment, method = "draws"
)$cl


assignment_nested <- mcclust.ext::minVI(
  psm = mcclust::comp.psm(res_nested$assignment),
  cls.draw = res_nested$assignment, method = "draws"
)$cl

print("Number of supernodes with coarsening of the likelihood:")
K_coars <- length(unique(assignment_coars))
print(K_coars)

print("Number of supernodes with nested MCMC:")
K_nested <- length(unique(assignment_nested))
print(K_nested)



print("Running GO overrepresentation analysis...")

# Reorder partition such that supernodes are indexed from larger to small.
tmp <- as.integer(names(sort(table(assignment_nested), decreasing = TRUE)))

supernode_nested <- factor(x = integer(p), levels = formatC(
  x = 1:K_nested, width = 2, format = "d", flag = "0"
))

for (k in 1:K_nested) supernode_nested[assignment_nested == tmp[k]] <- formatC(
  x = k, width = 2, format = "d", flag = "0"
)


tmp <- as.integer(names(sort(table(assignment_coars), decreasing = TRUE)))

supernode_coars <- factor(x = integer(p), levels = formatC(
  x = 1:K_coars, width = 2, format = "d", flag = "0"
))

for (k in 1:K_coars) supernode_coars[assignment_coars == tmp[k]] <- formatC(
  x = k, width = 2, format = "d", flag = "0"
)


library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

xx_both <- compareCluster(
  geneClusters = geneID ~ supernode + nested, data = data.frame(
    geneID = all_genes, supernode = c(supernode_nested, supernode_coars),
    nested = c(
      rep(TRUE, length(supernode_nested)), rep(FALSE, length(supernode_coars))
    )
  ), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
  universe = all_genes
)

xx_both@compareClusterResult$type <- "Coarsened likelihood"

xx_both@compareClusterResult$type[
  xx_both@compareClusterResult$nested == "TRUE"
] <- "Nested MCMC"

xx_both@compareClusterResult$supernode2 <- as.factor(
  as.integer(xx_both@compareClusterResult$supernode)
)

dotplot(
  object = xx_both, x = "supernode2", showCategory = 3L, font.size = 11,
  label_format = 100
) + facet_grid(
  rows = ~ type, scales = "free", space = "free"
) + xlab("Supernode") + ylab("GO term") + theme(
  panel.grid.major.x = element_blank(), legend.position = "bottom",
  legend.key.width = unit(0.8, "cm"),
  strip.background = element_rect(fill = "black"),
  strip.text = element_text(face = "bold", color = "white")
) + scale_colour_viridis_c(
  limit = c(0, 0.05), end = 0.98, direction = -1
) + labs(
  color =  substitute(paste("Adjusted ", italic("p"), "-value")),
  size = "Gene ratio"
) + scale_size_continuous(range = c(2, 5.5))

ggsave(filename = "GOenrichment.pdf", width = 11.5, height = 10)



print("Plotting the graph sphere...")
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



# Estimate the supergraph with the tessellation fixed.


# Function that runs MCMC and generates output for the supergraph only
run_MCMC_G_only <- function (
    X, assignment, n_iter, burnin = NULL
) {
  if (is.null(burnin)) burnin <- n_iter %/% 10L
  n <- nrow(X)
  n_centers <- max(assignment)
  G_MCMC <- array(data = NA_integer_, dim = c(n_iter, n_centers, n_centers))
  super_edge_mat <- matrix(data = 0L, nrow = n_centers, ncol = n_centers)
  U <- compute_U_unaugmented(X = X, assignment = assignment)
  
  # Supergraph: initialize at empty
  G <- matrix(0L, nrow = n_centers, ncol = n_centers)
  
  pb <- txtProgressBar(max = burnin + n_iter, style = 3)
  
  for (s in 1:(burnin + n_iter)) {
    G <- update_G_WWA(G, 0.5, 3, U, n)
    
    if (s > burnin) {
      G_MCMC[s - burnin, , ] <- G
      super_edge_mat <- super_edge_mat + G
    }
    
    setTxtProgressBar(pb, s)
  }
  
  close(pb)
  
  return(list(
    G_MCMC = G_MCMC, superedge_mat = super_edge_mat / n_iter
  ))
}



print("Running MCMC on the supergraph with the tessellation fixed...")
set.seed(1L)
res_G <- run_MCMC_G_only(X = X, assignment = assignment_nested, n_iter = 1e5L)
supergraph_adj <- res_G$superedge_mat > 0.99

print("Point estimate of the supergraph with 0.99 edge inclusion threshold:")
print(igraph::graph_from_adjacency_matrix(
  adjmatrix = supergraph_adj, mode = "undirected", diag = FALSE
))

p_k_sub <- sort(table(assignment_nested))[K_nested:1]
p_sub <- sum(p_k_sub)
reorder_supernodes <- as.integer(names(p_k_sub))

super_layout <- graphlayouts::layout_with_stress(
  g = igraph::graph_from_adjacency_matrix(
    adjmatrix = supergraph_adj[reorder_supernodes, reorder_supernodes],
    mode = "undirected", diag = FALSE
  ), tol = 1e-8, mds = TRUE
)

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
    15 * super_layout[k, ] + t(igraph::layout_as_tree(
      graph = G_list[[k]], circular = TRUE
    ))
  )
}

G_block <- igraph:::union.igraph(G_list, byname = TRUE)
tmp <- list(membership = assignment_nested[as.integer(igraph::V(G_block)$name)])
class(tmp) <- "communities"

cols <- RColorBrewer::brewer.pal(
  n = 6L, name = "Set1"
)[c(3L, 1L, 2L, 4L, 5L, 6L)]

mod_ind <- rep(1:6, mod_sizes)

pdf("gene_graph_sphere_nested.pdf", width = 8, height = 8)
par(mar = c(0, 0, 0, 0))

plot.igraph2(
  x = G_block, supergraph_adj = supergraph_adj, mark.groups = tmp,
  mark.col = scales::alpha("white", 0.5), mark.border = "darkgray",
  mark.expand = 6, vertex.size = 1.5,
  vertex.color = cols[mod_ind[as.integer(V(G_block)$name)]],
  vertex.label = "", edge.color = "black", layout = layout[, 2:1],
  xlim = c(-1.01, 1)
)

legend(
  x = 0.65, y = -0.5, legend = paste("Module", c(5, 6, 14, 27, 36, 39)),
  box.col = "white", pch = 21L, pt.bg = cols, cex = 1
)

dev.off()



print("Plotting graphical lasso estimate...")
pdf("gene_glasso.pdf", width = 16, height = 8)
par(mar = c(0, 0, 1, 0), mfrow = c(1L, 2L))

for (i in 1:2) {
  igraph::plot.igraph(
    x = igraph::graph_from_adjacency_matrix(adjmatrix = glasso::glasso(
      s = crossprod(X) / nrow(X), rho = 0.5568  # No. of edges matches `G_block`.
    )$wi != 0, mode = "undirected", diag = FALSE), vertex.size = 1.5,
    vertex.color = cols[mod_ind[as.integer(V(G_block)$name)]],
    vertex.label = "", edge.color = "black", layout = if (i == 1) layout[, 2:1],
    xlim = c(-1.01, 1), main = c("Fixed layout", "Automatic layout")[i]
  )
  
  if (i == 1) legend(
    x = 0.65, y = -0.5, legend = paste("Module", c(5, 6, 14, 27, 36, 39)),
    box.col = "white", pch = 21L, pt.bg = cols, cex = 1
  )
}

dev.off()



print(
  "Checking whether superedges appear more between supernodes corresponding to the same module..."
)

n_e <- sum(supergraph_adj) %/% 2L  # Number of edges (`n_e`) in the supergraph
n_e_possible_within <- 0L
n_e_within <- 0L


# Function that computes the module corresponding to the kth supernode
get_mod <- function (k) as.integer(names(sort(
  x = table(mod_ind[assignment_nested == k]), decreasing = TRUE
))[1])

for (k1 in 1:(K_nested - 1L)) {
  mod1 <- get_mod(k1)
  
  for (k2 in (k1 + 1L):K_nested) if (mod1 == get_mod(k2)) {
    n_e_possible_within <- n_e_possible_within + 1L
    if (supergraph_adj[k1, k2] == 1L) n_e_within <- n_e_within + 1L
  }
}


n_e_between <- n_e - n_e_within

n_e_possible_between <- (
  K_nested * (K_nested - 1L)
) %/% 2L - n_e_possible_within

print("Proportion of pairs of nodes with a superedge")
print("Within modules:")
print(n_e_within / n_e_possible_within)
print("Between modules:")
print(n_e_between / n_e_possible_between)



print("Inspecting gene interactions in the graph sphere...")
string_db <- STRINGdb::STRINGdb$new(version = "12.0", network_type = "full")

mapped <- string_db$map(
  my_data_frame = data.frame(gene = all_genes),
  my_data_frame_id_col_names = "gene", takeFirst = TRUE,
  removeUnmappedRows = FALSE
)

tmp <- string_db$get_interactions(mapped$STRING_id)
tmp$from <- all_genes[match(tmp$from, mapped$STRING_id)]
tmp$to <- all_genes[match(tmp$to, mapped$STRING_id)]
inter_adj <- matrix(nrow = p, ncol = p)

for (i in 1:p) {
  ind <- which(tmp$from == all_genes[i])
  
  for (j in 1:p) if (!isTRUE(inter_adj[i, j])) {
    inter_adj[i, j] <- any(tmp$to[ind] == all_genes[j])
  }
}

for (i in 1:p) for (j in 1:p) {
  inter_adj[i, j] <- inter_adj[i, j] | inter_adj[j, i]
}

prop_between <- matrix(data = NA_real_, nrow = K_nested, ncol = K_nested)

for (k1 in 1:K_nested) for (k2 in k1:K_nested) prop_between[k1, k2] <- mean(
  inter_adj[assignment_nested == k1, assignment_nested == k2]
)

prop_between_no_diag <- prop_between
diag(prop_between_no_diag) <- NA_real_

levels <- c(
  "Within supernode", "Between supernodes:\nWith superedge",
  "Between supernodes:\nWithout superedge"
)

pdf("gene_interactions.pdf", width = 5, height = 10)
par(las = 2, mar = c(10, 4, 4, 2) + 0.1)

boxplot(
  formula = prop ~ label, data = data.frame(prop = c(
    diag(prop_between), prop_between[supergraph_adj],
    prop_between_no_diag[!supergraph_adj]
  ), label = factor(x = rep(levels, c(
    K_nested, sum(supergraph_adj), sum(!supergraph_adj)
  )), levels = levels)), range = 0, xlab = "",
  ylab = "Proportion of gene pairs with an interaction", ylim = c(0, 0.25)
)

dev.off()

max(prop_between, na.rm = TRUE)
max(prop_between_no_diag, na.rm = TRUE)
