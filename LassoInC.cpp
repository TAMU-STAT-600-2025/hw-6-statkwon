#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda) {
  return ((a > 0) - (a < 0)) * std::max(abs(a) - lambda, 0.);
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde,
               const arma::colvec& beta, double lambda) {
  auto n = Xtilde.n_rows;
  auto res = Ytilde - (Xtilde * beta);
  double f_obj = arma::as_scalar(res.t() * res) / (2 * n) +
    lambda * arma::accu(arma::abs(beta));
  return f_obj;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde,
                                    const arma::colvec& Ytilde, double lambda,
                                    const arma::colvec& beta_start,
                                    double eps = 0.001) {
  auto n = Xtilde.n_rows;
  auto p = Xtilde.n_cols;
  
  arma::colvec beta_start_ = beta_start;
  if (beta_start_.is_empty()) {
    beta_start_ = arma::zeros(p);
  }
  
  arma::colvec r = Ytilde - Xtilde * beta_start_;
  arma::colvec beta = beta_start_;
  double error = 1000.;
  while (error > eps) {
    arma::colvec beta_old = beta;
    for (int j = 0; j < p; j++) {
      beta[j] = soft_c(beta_old[j] + arma::as_scalar(Xtilde.col(j).t() * r) / n,
                       lambda);
      r = r + Xtilde.col(j) * (beta_old[j] - beta[j]);
    }
    error = lasso_c(Xtilde, Ytilde, beta_old, lambda) -
      lasso_c(Xtilde, Ytilde, beta, lambda);
  }
  
  return beta;
}

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde,
                                     const arma::colvec& Ytilde,
                                     const arma::colvec& lambda_seq,
                                     double eps = 0.001) {
  auto p = Xtilde.n_cols;
  auto n_lambda = lambda_seq.size();
  
  arma::mat betas(p, n_lambda);
  arma::colvec beta_start = arma::zeros(p);
  
  for (int j = 0; j < n_lambda; j++) {
    betas.col(j) =
      fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq[j], beta_start, eps);
    beta_start = betas.col(j);
  }
  
  return betas;
}
