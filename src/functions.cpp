#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double negllk_norm(
                   NumericVector Betas,
                   NumericMatrix X,
                   NumericVector y,
                   NumericVector Snewdiag,
                   double lambda
                   ){
  int n = y.length();
  int p = Betas.length(); 

  double sum = 0.0;
  for(int i = 0; i < n; i ++){
    double Xb = 0.0;
    for(int j = 0; j < p; j ++){
      Xb = Xb + X(i,j) * Betas[j];
    }
    sum = sum + (y[i] - Xb) * (y[i] - Xb);
  }
  
  double pen = 0.0;
  for (int i = 0; i < p; i++) {
      pen = pen + Snewdiag[i] * Betas[i] * Betas[i];
  }
  
  
  return(sum + pen * lambda);
}