# graph-sphere

Repository with code used for the paper "Graph Sphere: From Nodes to Supernodes
in Graphical Models" by Willem van den Boom, Maria De Iorio, Alexandros Beskos
and Ajay Jasra ([arXiv:2310.11741](https://arxiv.org/abs/2310.11741))


## Description of files

* [`graph_sphere_MCMC.R`](graph_sphere_MCMC.R) implements the Markov chain
Monte Carlo (MCMC) for the graph sphere model.
[`graph_sphere_MCMC.R`](graph_sphere_MCMC.R) compiles the C++ code in the
folder [src](src/).

* [`gene.R`](gene.R) contains code for the application with gene expression
data. [`gene.R`](gene.R) loads [`graph_sphere_MCMC.R`](graph_sphere_MCMC.R).

* [`simulations.R`](simulations.R) contains ode for the simulation studies in
Appendix H of Supplementary Material. [`simulations.R`](simulations.R) loads
[`graph_sphere_MCMC.R`](graph_sphere_MCMC.R).

* The folder [src](src/) contains C++ scripts for MCMC and a Laplace
approximation for Bayesian Gaussian graphical models (GGMs) with the
*G*-Wishart distribution as prior on the precision matrix.
