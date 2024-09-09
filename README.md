# graph-of-graphs

Repository with code used for the paper "Graph of Graphs: From Nodes to
Supernodes in Graphical Models" by Maria De Iorio, Willem van den Boom,
Alexandros Beskos, Ajay Jasra and Andrea Cremaschi
([arXiv:2310.11741](https://arxiv.org/abs/2310.11741))


## Description of files

* [`graph_of_graphs_MCMC.R`](graph_of_graphs_MCMC.R) implements the Markov
chain Monte Carlo (MCMC) for the graph of graphs model.
[`graph_of_graphs_MCMC.R`](graph_of_graphs_MCMC.R) compiles the C++ code in the
folder [`src`](src/).

* [`gene.R`](gene.R) contains code for the application with gene expression
data. [`gene.R`](gene.R) loads
[`graph_of_graphs_MCMC.R`](graph_of_graphs_MCMC.R).

* The folder [`simul`](simul/) contains code for the simulation studies in
Supplementary Material. [`simul_lat_fac.R`](simul/simul_lat_fac.R) and
[`simul_mvt.R`](simul/simul_mvt.R) load
[`graph_of_graphs_MCMC.R`](graph_of_graphs_MCMC.R).

* The folder [`src`](src/) contains C++ scripts for MCMC and a Laplace
approximation for Bayesian Gaussian graphical models (GGMs) with the
*G*-Wishart distribution as prior on the precision matrix.

* The folder [`misc`](misc/) contains scripts for miscellaneous other plots.

* [`run_all.sh`](run_all.sh) is a shell script that provides a workflow for how
to run the code.
