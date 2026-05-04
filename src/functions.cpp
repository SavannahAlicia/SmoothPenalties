#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double negllk_norm(
                   NumericVector Betas,
                   NumericMatrix X,
                   NumericVector y,
                   NumericMatrix S,
                   double lambda
                   ){
  int n = y.length();
  int p = Betas.length() - 1; // sigma is last beta
  double sigma = Betas[(p)];
  double lhs = (n/2) * log(2 * M_PI * sigma* sigma);
  double sum = 0.0;
  for(int i = 0; i < n; i ++){
    double Xb = 0.0;
    for(int j = 0; j < p; j ++){
      Xb = Xb + X(i,j) * Betas[j];
    }
    sum = sum + (y[i] - Xb) * (y[i] - Xb);
  }
  double rhs = (1/2 * sigma * sigma) * sum;
  
  double pen = 0.0;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      pen = pen + Betas[i] * S(i, j) * Betas[j];
    }
  }
  
  pen *= 0.5 * lambda;
  
  return(lhs + rhs + pen);
}