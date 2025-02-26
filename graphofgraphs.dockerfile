FROM r-base:4.4.2

RUN R -e 'install.packages("remotes")'
RUN R -e 'library(remotes); install_version("abind")'
RUN R -e 'library(remotes); install_version("BH")'
RUN R -e 'library(remotes); install_version("CholWishart")'
RUN R -e 'library(remotes); install_version("latticeExtra")'
RUN R -e 'library(remotes); install_version("matrixStats")'
RUN R -e 'library(remotes); install_version("memoise")'
RUN R -e 'library(remotes); install_version("mvtnorm")'
RUN R -e 'library(remotes); install_version("RcppBlaze")'
RUN R -e 'library(remotes); install_version("withr")'

RUN R -e 'withr::with_makevars(new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"), code = Rcpp::sourceCpp(file = "src/ggm.cpp"))'






