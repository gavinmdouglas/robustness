// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP fast_matrix_multiply(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
	Eigen::MatrixXd C = A * B;

	return(Rcpp::wrap(C));
}