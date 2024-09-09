### Graph of graphs: Visualization for Proposition 4

# Code for the visualization of the lower bound from Proposition 4(i) in the
# paper "Graph of Graphs: From Nodes to Supernodes in Graphical Models" by Maria
# De Iorio, Willem van den Boom, Alexandros Beskos, Ajay Jasra and Andrea
# Cremaschi.


library(ggplot2)


p <- 10L
h_p <- function (t) 1 / p + (1 - 1 / p) * t

h_star <- function (t) {
  n <- length(t)
  res <- numeric(n)
  
  for (i in 1:n) {
    res[i] <- if (t[i] <= 0.5) t[i] else 0.5 * (1 + sqrt(2 * t[i] - 1))
  }
  
  return (res)
}

h_p_inv <- function (t) (p * t - 1) / (p - 1L)
h_star_inv <- function (t) ifelse(t <= 0.5, t, 0.5 * (1 + (2 * t - 1)^2))

n_grid <- 1e3L
t_grid <- seq(from = 0, to = 1, length.out = n_grid)
df <- data.frame(rho = numeric(n_grid^2), s2 = numeric(n_grid^2))
pb <- txtProgressBar(max = n_grid, style = 3L)

for (i in 1:n_grid) {
  df$rho[0:(n_grid - 1L) * n_grid + i] <- t_grid
  df$s2[(i - 1L) * n_grid + 1:n_grid] <- t_grid
  setTxtProgressBar(pb, i)
}

close(pb)
df$lb <- pmax(h_p(df$rho), h_star(h_p(df$s2)))

plot_obj <- ggplot() +
  geom_raster(data = df, aes(x = rho, y = s2, fill = lb)) +
  scale_fill_viridis_c(name = "Lower bound", limits = 0:1) +
  geom_function(
    fun = function (rho) h_p_inv(h_star_inv(h_p(rho))), col = "red",
    lty = "dashed", lwd = 1
  ) +
  coord_fixed() +
  scale_x_continuous(limits = 0:1, expand = c(0, 0)) +
  scale_y_continuous(limits = 0:1, expand = c(0, 0)) +
  theme(
    panel.border = element_rect(colour = "black", fill= NA, linewidth = 1),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  ylab(expression(italic(s)^2)) +
  xlab(expression(rho))

# We suppress warnings related to values not being plot since they are exactly
# at the plot limits.
suppressWarnings(ggsave(
  plot = plot_obj, filename = "proposition_plot.pdf", width = 5.5, height = 4
))
