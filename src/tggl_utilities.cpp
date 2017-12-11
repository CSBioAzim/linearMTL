#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat RcppShrink(arma::mat A, arma::mat ranges) {
/*
 * Applies shrinkage operator to matrix A.
 *
 * A is matrix of size (cumsum(inner group sizes)) by J.
 * Shrink all vectors A[group ids, j].
 */
  int V = ranges.n_rows;
  int v, k;
  arma::rowvec gnorm;
  for (v = 0; v < V; v++) {
    gnorm = arma::zeros<arma::rowvec>(A.n_cols);
    for (k = 0; k < A.n_cols; k++) {
      gnorm(k) = arma::accu(pow(A(arma::span(ranges(v, 0)-1, ranges(v, 1)-1), k), 2));
    }
    gnorm = arma::sqrt(gnorm);
    for (k = 0; k < A.n_cols; k++) {
      A(arma::span(ranges(v, 0)-1, ranges(v, 1)-1), k) = A(arma::span(ranges(v, 0)-1, ranges(v, 1)-1), k) / std::max(gnorm(k), 1.0);
    }
  }
  return A;
}
