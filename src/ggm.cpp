// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppBlaze)]]
// [[Rcpp::plugins(cpp17)]]

// OpenMP can be enabled by "// [[Rcpp::plugins(openmp)]]" but requires
// additional compiler flags on macOS, and the OpenMP code is commented out
// (see line 398) due to an OpenMP issue on Windows when using RcppBlaze.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Avoid OpenMP linker error on Windows:
#define BLAZE_USE_SHARED_MEMORY_PARALLELIZATION 0
#include <RcppBlaze.h>

// The distributions in `<random>` are not portable. That is, they do not
// yield the same random numbers on different machines. Therefore, we use the
// distributions from Boost, which are protable.
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>



// C++'s RNGs are not very fast. Therefore, I use the RNG from
// https://gist.github.com/martinus/c43d99ad0008e11fcdbf06982e25f464:
// extremely fast random number generator that also produces very high quality random.
// see PractRand: http://pracrand.sourceforge.net/PractRand.txt
class sfc64 {
  public:
    using result_type = uint64_t;

    static constexpr uint64_t(min)() { return 0; }
    static constexpr uint64_t(max)() { return UINT64_C(-1); }

    sfc64() : sfc64(std::random_device{}()) {}

    explicit sfc64(uint64_t seed) : m_a(seed), m_b(seed), m_c(seed), m_counter(1) {
        for (int i = 0; i < 12; ++i) {
            operator()();
        }
    }

    uint64_t operator()() noexcept {
        auto const tmp = m_a + m_b + m_counter++;
        m_a = m_b ^ (m_b >> right_shift);
        m_b = m_c + (m_c << left_shift);
        m_c = rotl(m_c, rotation) + tmp;
        return tmp;
    }

  private:
    template <typename T> T rotl(T const x, int k) { return (x << k) | (x >> (8 * sizeof(T) - k)); }

    static constexpr int rotation = 24;
    static constexpr int right_shift = 11;
    static constexpr int left_shift = 3;
    uint64_t m_a;
    uint64_t m_b;
    uint64_t m_c;
    uint64_t m_counter;
};


template <typename T>
blaze::DynamicMatrix<T> R2Blaze(Rcpp::NumericMatrix& mat_R) {
    // Convert an R matrix to a Blaze matrix.
    int n_row = mat_R.nrow(), n_col = mat_R.ncol();
    blaze::DynamicMatrix<T> mat_Blaze(n_row, n_col);

    for (int i = 0; i < n_row; i++) for (int j = 0; j < n_col; j++)
        mat_Blaze(i, j) = mat_R(i, j);

    return mat_Blaze;
}


template <typename T>
Rcpp::NumericMatrix Blaze2R(blaze::DynamicMatrix<T>& mat_Blaze) {
    // Convert a Blaze matrix to an R matrix.
    int n_row = mat_Blaze.rows(), n_col = mat_Blaze.columns();
    Rcpp::NumericMatrix mat_R(n_row, n_col);

    for (int i = 0; i < n_row; i++) for (int j = 0; j < n_col; j++)
        mat_R(i, j) = mat_Blaze(i, j);

    return mat_R;
}


blaze::DynamicMatrix<double> inv_pos_def(blaze::DynamicMatrix<double> mat) {
    blaze::invert<blaze::byLLH>(mat);
    return mat;
}


blaze::UpperMatrix<blaze::DynamicMatrix<double> > rwish_identity_chol(
    int p, double df, sfc64& rng
) {
    blaze::UpperMatrix<blaze::DynamicMatrix<double> > Phi(p);
    boost::random::normal_distribution<> rnorm(0.0, 1.0);
    df += p - 1;
    
    // Generate the upper-triangular Cholesky decompositon of a standard
    // Wishart random variable.
    for (int i = 0; i < p; i++) {
        boost::random::chi_squared_distribution<> rchisq(df - i);
        Phi(i, i) = std::sqrt(rchisq(rng));
        for (int j = i + 1; j < p; j++) Phi(i, j) = rnorm(rng);
    }
    
    return Phi;
}


blaze::DynamicMatrix<double> rwish_identity(int p, double df, sfc64& rng) {
    /*
    Sample a `p` by `p` matrix from a Wishart distribution with `df` degrees of
    freedom with an identity rate matrix.
    */

    blaze::UpperMatrix<blaze::DynamicMatrix<double> >
        Phi = rwish_identity_chol(p, df, rng);
    
    return declsym(trans(Phi) * Phi);
}


blaze::DynamicMatrix<double> rwish(
    double df, blaze::DynamicMatrix<double>& rate, sfc64& rng
) {
    blaze::UpperMatrix<blaze::DynamicMatrix<double> >
        Phi_identity = rwish_identity_chol(rate.rows(), df, rng);

    blaze::LowerMatrix<blaze::DynamicMatrix<double> > chol;
    llh(rate, chol);
    blaze::invert(chol);
    blaze::DynamicMatrix<double> Phi = Phi_identity * chol;
    return declsym(trans(Phi) * Phi);
}


blaze::DynamicMatrix<double> rgwish_body(
    blaze::DynamicMatrix<int>& adj,
    blaze::SymmetricMatrix<blaze::DynamicMatrix<double> > W
) {
    int p = adj.rows();
    if (sum(adj) == p * (p - 1)) return W;  // The graph is complete.

    // This function follows Section 2.4 in Lenkoski (2013, arXiv:1304.1350v1).
    blaze::SymmetricMatrix<blaze::DynamicMatrix<double> >
        Sigma = inv_pos_def(W);

    W = Sigma;  // Step 1

    // Inspired by C++ code from the R package BDgraph.
    // Avoid recomputing the neighbors for each iteration:
    std::vector<std::vector<double> > neighbors(p);
    std::vector<blaze::DynamicVector<double> > Sigma_N(p);

    for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) {
        if (adj(i, j)) neighbors[i].push_back(j);
        Sigma_N[i] = elements(column(Sigma, i), neighbors[i]);
    }

    blaze::SymmetricMatrix<blaze::DynamicMatrix<double> > W_previous(p);
    blaze::DynamicMatrix<double, blaze::columnMajor> W_N;
    blaze::DynamicMatrix<double, blaze::rowMajor> W_NN;
    blaze::DynamicVector<double> W_beta_hat(p), beta_star;
    
    for (int i = 0; i < 1000; i++) {
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
            // The next line is not needed as Blaze enforces symmetry of `W`.
            // row(W, j) = trans(W_beta_hat);
            W(j, j) = W_jj;
        }

        // 1e-8 is consistent with BDgraph::rgwish.
        if (blaze::mean(blaze::abs(W - W_previous)) < 1e-5)
            return inv_pos_def(W);
    }

    std::cout << "`rgwish` failed to converge." << std::endl;
    return inv_pos_def(W);
}


blaze::DynamicMatrix<double> rgwish_identity(
    blaze::DynamicMatrix<int>& adj, double df, sfc64& rng
) {
    /*
    Sample from the G-Wishart.
    
    This function assumes an identity scale matrix.
    */
    return rgwish_body(adj, rwish_identity(adj.rows(), df, rng));
}


blaze::DynamicMatrix<double> rgwish(
    blaze::DynamicMatrix<int>& adj, double df,
    blaze::DynamicMatrix<double>& rate, sfc64& rng
) {
    /*
    Sample from the G-Wishart using the Lenkoski method.
    
    `rate` is the inverse of the scale matrix of the Wishart distribution.
    This function follows Section 2.4 in Lenkoski (2013, arXiv:1304.1350v1).
    */
    return rgwish_body(adj, rwish(df, rate, rng));
}


double proposal_G_es(int p, int n_e_tilde, int n_e) {
    // Proposal transition probability from `G` to `G_tilde` based on edge
    // counts
    int max_e = p * (p - 1) / 2;

    if (n_e == 0 or n_e == max_e) {
        return 1.0 / max_e;
    } else if (n_e > n_e_tilde) {
        return 0.5 / n_e;
    } else {
        return 0.5 / (max_e - n_e);
    }
}


template <typename T>
blaze::DynamicMatrix<T> permute_mat(
    blaze::DynamicMatrix<T>& mat, std::vector<int>& perm
) {
    std::vector<int> from(0), to(0);

    for (int i = 0; i < perm.size(); i++) if (perm[i] != i) {
        from.push_back(perm[i]);
        to.push_back(i);
    }

    blaze::DynamicMatrix<T> perm_rows = rows(mat, from);
    columns(perm_rows, to) = columns(perm_rows, from);
    blaze::DynamicMatrix<T> mat_perm(mat);
    rows(mat_perm, to) = perm_rows;
    columns(mat_perm, to) = trans(perm_rows);
    return mat_perm;
}


double log_N_tilde(
    blaze::LowerMatrix<blaze::DynamicMatrix<double> >& Phi,
    double rate_perm_11, double rate_perm_21
) {
    /*
    The log of the function N from
    Cheng & Lengkoski (2012, page 2314, doi:10.1214/12-EJS746)

    `rate_perm_11` and `rate_perm_21` contain the element in the last row and
    column and the element just above it in the rate matrix, respectively.
    */
    int p = Phi.rows();

    return std::log(Phi(p - 2, p - 2)) + 0.5*(
        -std::log(rate_perm_11) + rate_perm_11*std::pow(
            -sum(
                submatrix(Phi, p - 2, 0, 1, p - 2)
                    % submatrix(Phi, p - 1, 0, 1, p - 2)
            // The plus in the next line is a minus in Cheng & Lengkoski
            // (2012). I believe that the minus is a typo in the article.
            )/Phi(p - 2, p - 2) + Phi(p - 2, p - 2)*rate_perm_21/rate_perm_11,
            2
        )
    );
}


double log_norm_ratio_Letac(
    blaze::DynamicMatrix<int>& adj, int i, int j, double df_0
) {
    /*
    Log of the approximation of the ratio of normalizing constants of the
    G-Wishart prior distributions from Letac et al. (2018, arXiv:1706.04416v2).

    The ratio is evaluated at the graph given by adjacency matrix `adj` with
    edge (`i`, `j`) absent (numerator) divided by the same graph with
    (`i`, `j`) present (denominator).

    `df_0` is the degrees of freedom of the G-Wishart prior distribution.
    */
    // `n_paths` is the number of paths of length 2 that connect nodes `i` and
    // `j`.
    int p = adj.rows(), n_paths = 0;
    for (int l = 0; l < p; l++) if (adj(i, l) and adj(j, l)) n_paths++;

    return std::log(0.5) - 0.5 * std::log(M_PI)
        + std::lgamma(0.5 * (df_0 + n_paths))
        - std::lgamma(0.5 * (df_0 + n_paths + 1.0));
}


std::tuple<std::vector<int>, std::vector<int>> permute_e_last(
    int i, int j, int p
) {
    /*
    Permute the nodes such that edge `(i, j)` becomes edges (p - 1, p).

    This function returns the permutation and inverse permutation.
    */
    std::vector<int> perm(p), perm_inv(p);
    std::iota(perm.begin(), perm.end(), 0);

    // Permute the nodes involved in `e`.
    if (i != p - 2) {
        perm[i] = p - 2;

        if (j == p - 2) {
            perm[p - 2] = p - 1;
            perm[p - 1] = i;
        } else {
            perm[p - 2] = i;
            perm[p - 1] = j;
            perm[j] = p - 1;
        }
    }

    for (int l = 0; l < p; l++) perm_inv[perm[l]] = l;    
    return std::make_tuple(perm, perm_inv);
}


double log_balancing_function(double log_t) {
    /*
    Compute the log of the balancing function t/(1+t).

    This function equals -log1pexp(-`log_t`). We use Equation 10 from
    https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
    to compute it.
    */
    // return 0.5 * log_t;
    // if (burnin) return 0.5 * log_t;
    if (log_t > 37.0) return -std::exp(-log_t);
    if (log_t > -18.0) return -std::log1p(std::exp(-log_t));
    if (log_t > -33.3) return log_t - std::exp(log_t);
    return log_t;
}


std::tuple<
    blaze::DynamicMatrix<double>, blaze::DynamicMatrix<double>
> locally_balanced_proposal(
    blaze::DynamicMatrix<double>& K, blaze::DynamicMatrix<int>& adj, int n_e,
    blaze::DynamicMatrix<double>& edge_prob_mat, double df_0,
    blaze::DynamicMatrix<double>& rate
) {
    /*
    Compute the locally balanced proposal from
    Zanella (2019, doi:10.1080/01621459.2019.1585255)
    */
    int p = adj.rows();
    blaze::DynamicMatrix<double> log_Q(p, p, -INFINITY);
    // The matrix `Phi` is specified here to avoid repeated expensive memory
    // (de)allocation.
    blaze::LowerMatrix<blaze::DynamicMatrix<double> > Phi;

    double log_q_add = std::log(proposal_G_es(p, n_e + 1, n_e)),
        log_q_rm = std::log(proposal_G_es(p, n_e - 1, n_e));

    // Fused triangular loop based on 
    // https://stackoverflow.com/a/33836073/5216563 as OpenMP does not support
    // collapsing triangular loops.
    // Commented out to avoid OpenMP linker error with RcppBlaze on Windows.
    // #pragma omp parallel for schedule(static) private(Phi)
    for (int e_id = 0; e_id < p * (p - 1) / 2; e_id++) {
        int i = e_id / p, j = e_id % p;

        if (i >= j) {
            i = p - i - 2;
            j = p - j - 1;
        }

        auto [perm, perm_inv] = permute_e_last(i, j, p);
        blaze::DynamicMatrix<double> K_perm = permute_mat(K, perm_inv);
        llh(K_perm, Phi);
        int exponent;

        if (adj(i, j)) {
            exponent = -1;
            log_Q(i, j) = log_q_rm;
        } else {
            exponent = 1;
            log_Q(i, j) = log_q_add;
        }

        log_Q(i, j) += log_balancing_function(exponent * (
            std::log(edge_prob_mat(i, j)) - std::log1p(-edge_prob_mat(i, j))
                + log_N_tilde(
                    Phi, rate(perm_inv[p - 1], perm_inv[p - 1]),
                    rate(perm_inv[p - 2], perm_inv[p - 1])
                ) + log_norm_ratio_Letac(adj, i, j, df_0)
        ));
    }

    log_Q -= max(log_Q);
    blaze::DynamicMatrix<double> Q = exp(log_Q);
    double sum_Q = sum(Q);
    return std::make_tuple(Q / sum_Q, log_Q - std::log(sum_Q));
}


void update_K_from_Phi(
    int p, int j, std::vector<int>& perm, blaze::DynamicMatrix<double>& K,
    blaze::LowerMatrix<blaze::DynamicMatrix<double> >& Phi
) {
    /*
    Update the precision matrix `K` in place according to a new Phi.

    The matrix `Phi` has been permuted such that indices `i` and `j` become
    indices p-1 and p. Matrix `K` is unpermuted.
    */
    blaze::DynamicVector<double> K_vec = Phi * trans(row(Phi, p - 1));
    K_vec = elements(K_vec, perm);  // Undo the permutation.
    column(K, j) = K_vec;
    row(K, j) = trans(K_vec);
}


void update_single_edge(
    blaze::DynamicMatrix<double>& K, blaze::DynamicMatrix<int>& adj, int& n_e,
    blaze::DynamicMatrix<double>& edge_prob_mat, double df, double df_0,
    blaze::DynamicMatrix<double>& rate, sfc64& rng,
    blaze::LowerMatrix<blaze::DynamicMatrix<double> >& Phi,
    bool approx = false, bool delayed_accept = true, bool loc_bal = true
) {
    /*
    MCMC step that attempts to update a single edge

    This function modifies `K`, `adj` and `n_e` in place.

    The MCMC is similar to the one described in
    Cheng & Lenkoski (2012, doi:10.1214/12-EJS746).

    The matrix `Phi` is passed to avoid potentially expensive memory
    (de)allocation.

    `approx` indicates whether the target is approximated through the
    approximation of the ratio of normalizing constants of the G-Wishart prior
    distributions from Letac et al. (2018, arXiv:1706.04416v2).

    `delayed_accept` indicates whether the delayed acceptance MCMC in
    Algorithm 1 of Christen & Fox (2005, doi:10.1198/106186005X76983) is used.
    The surrogate posterior considered derives from the approximation of the
    ratio of normalizing constants of the G-Wishart prior distributions from
    Letac et al. (2018, arXiv:1706.04416v2).

    `approx` and `delayed_accept` cannot be true simultaneously.

    `loc_bal` indicates whether to use the locally balanced proposal from
    Zanella (2019, doi:10.1080/01621459.2019.1585255).
    */
    if (approx and delayed_accept) throw std::runtime_error(
        "`approx` and `delayed_accept` cannot be true simultaneously."
    );

    // (i, j) is the edge to be updated.
    int i, j, p = adj.rows(), max_e = p * (p - 1) / 2;
    bool accept = true, add;  // Whether an edge is being added
    // The steps and notation follow Algorithm 1 of Christen & Fox (2005).
    double log_q_y_x, log_q_x_y;  // Proposal transition probabilities
    boost::random::uniform_01<double> runif;
    
    // Sample from the proposal
    if (loc_bal) {
        // Compute the locally balanced proposal.
        auto [Q, log_Q] = locally_balanced_proposal(
            K, adj, n_e, edge_prob_mat, df_0, rate
        );

        // Sample an edge from the proposal using the inverse CDF method.
        i = 0;
        j = 0;
        double tmp_sum = 0.0, U = runif(rng);

        do {
            if (j == p - 1) j = ++i;
            tmp_sum += Q(i, ++j);
        } while (tmp_sum < U);

        add = not adj(i, j);
        log_q_y_x = log_Q(i, j);
    } else {
        // Decide whether to propose an edge addition or removal.
        if (n_e == 0) {
            add = true;
        } else if (n_e == max_e) {
            add = false;
        } else {
            add = runif(rng) < 0.5;
        }

        // Pick edge to add or remove uniformly at random.
        int row_sum, row_sum_cum = 0;
        i = -1;

        if (add) {
            boost::random::uniform_int_distribution<int>
                r_edge(1, max_e - n_e);

            int e_id = r_edge(rng);

            do {
                i++;

                row_sum
                    = p - i - 1 - sum(submatrix(adj, i, i + 1, 1, p - i - 1));
                
                row_sum_cum += row_sum;
            } while (row_sum_cum < e_id);

            e_id -= row_sum_cum - row_sum;
            j = i;
            row_sum = 0;
            while (row_sum < e_id) row_sum += 1 - adj(i, ++j);
        } else {
            boost::random::uniform_int_distribution<int> r_edge(1, n_e);
            int e_id = r_edge(rng);

            do {
                i++;
                row_sum = sum(submatrix(adj, i, i + 1, 1, p - i - 1));
                row_sum_cum += row_sum;
            } while (row_sum_cum < e_id);

            e_id -= row_sum_cum - row_sum;
            j = i;
            row_sum = 0;
            while (row_sum < e_id) row_sum += adj(i, ++j);
        }

        int n_e_tilde = n_e + 2*add - 1;
        log_q_y_x = std::log(proposal_G_es(p, n_e_tilde, n_e));
        log_q_x_y = std::log(proposal_G_es(p, n_e, n_e_tilde));
    }

    int exponent = 2*add - 1, n_e_tilde = n_e + exponent;
    auto [perm, perm_inv] = permute_e_last(i, j, p);
    blaze::DynamicMatrix<double> K_perm = permute_mat(K, perm_inv);
    llh(K_perm, Phi);

    double log_target_ratio_approx, log_g_x_y,
        rate_pp = rate(perm_inv[p - 1], perm_inv[p - 1]),
        rate_1p = rate(perm_inv[p - 2], perm_inv[p - 1]),
        log_N_tilde_post = log_N_tilde(Phi, rate_pp, rate_1p),
        log_prior_ratio = exponent * (
            std::log(edge_prob_mat(i, j)) - std::log1p(-edge_prob_mat(i, j))
        );

    double Phi_12_cur = Phi(p - 1, p - 2);
    boost::random::chi_squared_distribution<> rchisq(df);
    Phi(p - 1, p - 1) = std::sqrt(rchisq(rng) / rate_pp);
    boost::random::normal_distribution<> rnorm(0.0, 1.0);

    if (loc_bal) {
        // Compute the reverse transition probability `log_q_x_y`.

        // Update `Phi(p - 1, p - 2)` according to the proposal.
        if (add) {  // The graph contains (`i`, `j`) in the proposal.
            Phi(p - 1, p - 2) = rnorm(rng)/std::sqrt(rate_pp)
                - Phi(p - 2, p - 2)*rate_1p/rate_pp;
        } else {  // The graph does not contain (`i`, `j`) in the proposal.
            Phi(p - 1, p - 2) = -sum(
                submatrix(Phi, p - 2, 0, 1, p - 2)
                    % submatrix(Phi, p - 1, 0, 1, p - 2)
            ) / Phi(p - 2, p - 2);
        }

        // Make `adj` equal to the adjacency matrix of the proposed graph.
        adj(i, j) = add;
        adj(j, i) = add;
        update_K_from_Phi(p, j, perm, K, Phi);

        auto [Q_tilde, log_Q_tilde] = locally_balanced_proposal(
            K, adj, n_e_tilde, edge_prob_mat, df_0, rate
        );

        log_q_x_y = log_Q_tilde(i, j);
    }

    if (delayed_accept or approx) {
        log_target_ratio_approx = log_prior_ratio + exponent*(
            log_N_tilde_post + log_norm_ratio_Letac(adj, i, j, df_0)
        );

        log_g_x_y
            = std::min(0.0, log_q_x_y - log_q_y_x + log_target_ratio_approx);

        // Step 2
        accept = std::log(runif(rng)) < log_g_x_y;
    }

    if (accept and not approx) {
        double log_q_star_y_x ,log_q_star_x_y;

        if (delayed_accept) {
            log_q_star_y_x = log_g_x_y + log_q_y_x;

            log_q_star_x_y
                = std::min(log_q_x_y, log_q_y_x - log_target_ratio_approx);
        } else {
            log_q_star_y_x = log_q_y_x;
            log_q_star_x_y = log_q_x_y;
        }

        // Step 3
        // Exchange algorithm to avoid evaluation of normalization constants
        blaze::DynamicMatrix<int> adj_perm = permute_mat(adj, perm_inv);
        adj_perm(p - 2, p - 1) = add;
        adj_perm(p - 1, p - 2) = add;

        blaze::DynamicMatrix<double>
            K_0_tilde = rgwish_identity(adj_perm, df_0, rng);

        blaze::LowerMatrix<blaze::DynamicMatrix<double> > Phi_0_tilde;
        llh(K_0_tilde, Phi_0_tilde);

        // The log of the function N from
        // Cheng & Lengkoski (2012, page 2314, doi:10.1214/12-EJS746) with
        // identity rate matrix
        double log_N_tilde_prior
            = std::log(Phi_0_tilde(p - 2, p - 2)) + 0.5*std::pow(
                sum(
                    submatrix(Phi_0_tilde, p - 2, 0, 1, p - 2)
                        % submatrix(Phi_0_tilde, p - 1, 0, 1, p - 2)
                ) / Phi_0_tilde(p - 2, p - 2),
                2
            );

        double log_target_ratio = log_prior_ratio + exponent*(
            log_N_tilde_post - log_N_tilde_prior
        );

        accept = std::log(runif(rng))
            < log_q_star_x_y - log_q_star_y_x + log_target_ratio;
    }

    if (accept) {  // Update the graph.
        adj(i, j) = add;
        adj(j, i) = add;
        n_e = n_e_tilde;
    } else if (loc_bal) {  // Revert any update in `adj` and `Phi`.
        adj(i, j) = not add;
        adj(j, i) = not add;
        Phi(p - 1, p - 2) = Phi_12_cur;
    }

    if (not loc_bal) {
        // Update `Phi(p - 1, p - 2)`.
        if (adj(i, j)) {  // The graph contains (`i`, `j`) after updating.
            Phi(p - 1, p - 2) = rnorm(rng)/std::sqrt(rate_pp)
                - Phi(p - 2, p - 2)*rate_1p/rate_pp;
        } else {  // The graph does not contain (`i`, `j`) after updating.
            Phi(p - 1, p - 2) = -sum(
                submatrix(Phi, p - 2, 0, 1, p - 2)
                    % submatrix(Phi, p - 1, 0, 1, p - 2)
            ) / Phi(p - 2, p - 2);
        }
    }

    update_K_from_Phi(p, j, perm, K, Phi);  // Not required if `loc_bal and accept`?
}


// [[Rcpp::export]]
Rcpp::NumericMatrix update_G_Rcpp(
    Rcpp::NumericMatrix adj, Rcpp::NumericMatrix edge_prob_mat, double df,
    double df_0, Rcpp::NumericMatrix rate, int n_edge, int seed,
    bool approx = false, bool delayed_accept = true, bool loc_bal = true
) {
    /*
    `p` is the number of nodes.
    `adj_in` is the adjacency matrix of the graph.
    `n_edge` is the number of single edge MCMC steps that are attempted.

    `adj_in` is modified in place.

    Updating an `igraph_t` one edge at a time comes with a notable
    computational cost if `approx = true`. We therefore work with the
    adjacency matrix.
    */
    if (adj.nrow() == 1) return adj;  // No update if there is only one node
    sfc64 rng(seed);
    blaze::DynamicMatrix<int> adj_Blaze = R2Blaze<int>(adj);

    blaze::DynamicMatrix<double> rate_Blaze = R2Blaze<double>(rate),
        edge_prob_mat_Blaze = R2Blaze<double>(edge_prob_mat);

    int n_e = sum(adj_Blaze) / 2;
    blaze::DynamicMatrix<double> K = rgwish(adj_Blaze, df, rate_Blaze, rng);
    blaze::LowerMatrix<blaze::DynamicMatrix<double> > Phi(K.rows());

    for (int i = 0; i < n_edge; i++) update_single_edge(
        K, adj_Blaze, n_e, edge_prob_mat_Blaze, df, df_0, rate_Blaze, rng, Phi,
        approx, delayed_accept, loc_bal
    );

    return Blaze2R(adj_Blaze);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix rgwish_Rcpp(
    Rcpp::NumericMatrix adj, double df, Rcpp::NumericMatrix rate, int seed
) {
    // Function that enables the use of `rgwish()` from R.
    sfc64 rng(seed);
    blaze::DynamicMatrix<int> adj_Blaze = R2Blaze<int>(adj);

    blaze::DynamicMatrix<double> rate_Blaze = R2Blaze<double>(rate),
        K = rgwish(adj_Blaze, df, rate_Blaze, rng);

    return Blaze2R(K);
}
