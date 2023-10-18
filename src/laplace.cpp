#include <cmath>
#include <iostream>
#include <vector>

// This code was tested using Blaze version 3.8.0
// (https://bitbucket.org/blaze-lib/blaze/src/master/) with the fix from this
// pull request: https://bitbucket.org/blaze-lib/blaze/pull-requests/46.
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Column.h>
#include <blaze/math/Columns.h>
#include <blaze/math/Row.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Elements.h>
#include <blaze/math/Subvector.h>
#include <blaze/math/Band.h>


#include <Rcpp.h>
// [[Rcpp::plugins(cpp17)]]



blaze::DynamicMatrix<double> inv_pos_def(blaze::DynamicMatrix<double> mat) {
    blaze::invert<blaze::byLLH>(mat);
    return mat;
}


template <class T>
double log_det(T& mat) {
    // Log of the matrix determinant of `mat`
    blaze::DynamicMatrix<double> L;
    llh(mat, L);
    return 2.0 * sum(log(diagonal(L)));
}


blaze::DynamicMatrix<double> gwish_mode_inv(
    blaze::DynamicMatrix<int>& adj, double df,
    blaze::DynamicMatrix<double>& rate
) {
    /*
    Find the inverse of the mode of a G-Wishart distribution.
    
    `adj` is the adjacency matrix of the graph to which the precision is
    constrained.

    `df` is the degrees of freedom of the distribution.

    `rate` is the rate or inverse scale matrix of the distribution and must be
    symmetric positive definite.
    
    The notation in this function follows Section 2.4 in
    Lenkoski (2013, arXiv:1304.1350v1).
    
    The optimization procedure is presented in Algorithm 17.1 of
    the Elements of Statistical Learning by Hastie et al.
    
    Compare Equation 7 from
    https://proceedings.neurips.cc/paper/2009/hash/a1519de5b5d44b31a01de013b9b51a80-Abstract.html
    with Equation 17.11 of the Elements of Statistical Learning by
    Hastie et al. to understand why we set `Sigma = rate / (df - 2.0)`.
    */
    int p = rate.rows();
    blaze::DynamicMatrix<double> Sigma(rate / (df - 2.0)), W(Sigma);  // Step 1

    // Inspired by C++ code from the R package BDgraph
    // Avoid recomputing the neighbors for each iteration:
    std::vector<std::vector<double> > neighbors(p);
    std::vector<blaze::DynamicVector<double> > Sigma_N(p);

    for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) {
        if (adj(i, j)) neighbors[i].push_back(j);
        Sigma_N[i] = elements(column(Sigma, i), neighbors[i]);
    }

    blaze::DynamicMatrix<double> W_previous(p, p);
    blaze::DynamicMatrix<double, blaze::columnMajor> W_N;
    blaze::DynamicMatrix<double, blaze::rowMajor> W_NN;
    blaze::DynamicVector<double> W_beta_hat(p), beta_star;
    
    for (int i = 0; i < 10000; i++) {
        W_previous = W;
        
        for (int j = 0; j < p; j++) {
            if (neighbors[j].size() == 0) {
                // This only happens if the graph is not connected.
                W_beta_hat = 0.0;
            } else if (neighbors[j].size() == p - 1) {
                subvector(W_beta_hat, 0, j) = subvector(Sigma_N[j], 0, j);

                subvector(W_beta_hat, j + 1, p - j - 1)
                    = subvector(Sigma_N[j], j, p - j - 1);
            } else {
                W_N = columns(W, neighbors[j]);
                W_NN = rows(W_N, neighbors[j]);
                solve(declsym(W_NN), beta_star, Sigma_N[j]);
                W_beta_hat = W_N * beta_star;
            }

            double W_jj = W(j, j);
            column(W, j) = W_beta_hat;
            row(W, j) = trans(W_beta_hat);
            W(j, j) = W_jj;
        }

        // 1e-8 is consistent with BDgraph
        if (blaze::mean(blaze::abs(W - W_previous)) < 1e-8) return W;
    }

    std::cout << "`gwish_mode_inv` failed to converge." << std::endl;
    return W;
}


double log_gwish_norm_laplace(
    blaze::DynamicMatrix<int>& adj, double df,
    blaze::DynamicMatrix<double>& rate
) {
    /*
    Log of Laplace approximation of G-Wishart normalizing constant
    
    Log of the Laplace approximation of the normalizing constant of the
    G-Wishart distribution outlined by
    Lenkoski and Dobra (2011, doi:10.1198/jcgs.2010.08181)

    Only the diagonal of the Hessian matrix is used here, which is faster than
    a full Hessian but less accurate per Moghaddam et al. (2009,
    https://papers.nips.cc/paper/2009/hash/a1519de5b5d44b31a01de013b9b51a80-Abstract.html).
    */
    int p = rate.rows(), n_e = sum(adj) / 2, i = 0, j = 0;

    blaze::DynamicMatrix<double> K_inv = gwish_mode_inv(adj, df, rate),
        K = inv_pos_def(K_inv);

    double log_det_H = 2.0*sum(log(diagonal(K_inv))) + n_e*std::log(2.0),
        h = -0.5 * (trace(K * rate) - (df - 2.0)*log_det(K));

    for (int e_id = 0; e_id < n_e; e_id++) {
        while (adj(i, ++j) == 0) if (j == p - 1) j = ++i;

        log_det_H +=
            std::log(K_inv(i, i)*K_inv(j, j) + std::pow(K_inv(i, j), 2));
        
        if (j == p - 1) j = ++i;
    }

    // The sign of the Hessian `-0.5 * (df - 2.0) * H` is flipped compared to
    // Lenkoski and Dobra (2011, Section 4). I think that this is correct as
    // |Hessian| can be negative while |-Hessian| cannot.
    double ret = h + 0.5*(p + n_e)*std::log(2.0 * M_PI)
        - 0.5*((p + n_e)*std::log(0.5 * (df - 2.0)) + log_det_H);

    return ret;
}


blaze::DynamicMatrix<double> R2Blaze(Rcpp::NumericMatrix& mat_R) {
    // Convert an R matrix to a Blaze matrix.
    int n_row = mat_R.nrow(), n_col = mat_R.ncol();
    blaze::DynamicMatrix<double> mat_Blaze(n_row, n_col);

    for (int i = 0; i < n_row; i++) for (int j = 0; j < n_col; j++)
        mat_Blaze(i, j) = mat_R(i, j);

    return mat_Blaze;
}


blaze::DynamicMatrix<int> R2Blaze(Rcpp::IntegerMatrix& mat_R) {
    // Convert an R matrix to a Blaze matrix.
    int n_row = mat_R.nrow(), n_col = mat_R.ncol();
    blaze::DynamicMatrix<int> mat_Blaze(n_row, n_col);

    for (int i = 0; i < n_row; i++) for (int j = 0; j < n_col; j++)
        mat_Blaze(i, j) = mat_R(i, j);

    return mat_Blaze;
}


// [[Rcpp::export]]
double log_gwish_norm_laplace_Rcpp(
    Rcpp::IntegerMatrix adj, double df, Rcpp::NumericMatrix rate
) {
    blaze::DynamicMatrix<int> adj_Blaze = R2Blaze(adj);
    blaze::DynamicMatrix<double> rate_Blaze = R2Blaze(rate);
    return log_gwish_norm_laplace(adj_Blaze, df, rate_Blaze);
}
