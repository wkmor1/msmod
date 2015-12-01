#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the inner products of a predictor and coefficient matrix
//'
//' Given a matrix of predictors, X, and matrix of coefficients, B, calculate
//' the inner products, X[i, ] %*%  B[j, ].
//'
//' @param X numeric matrix. Predictors.
//' @param B numeric matrix. Coefficients. 
//' @export
// [[Rcpp::export]]

NumericMatrix inprod_mat(NumericMatrix X, NumericMatrix B) {
  
  int nX = X.nrow();
  int nB = B.nrow();
  int nK = X.ncol();
  
  NumericMatrix theta(nX, nB);
  
  for (int i = 0; i < nX; i++) {
    for (int j = 0; j < nB; j++) {
      theta(i, j) = 0;
      for (int k = 0; k < nK; k++) {
        theta(i, j) += X(i, k) * B(j, k);
      }
    }
  }

  return(theta);

}
