#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double negllk_norm(
                   NumericVector Betas,
                   NumericMatrix X,
                   NumericVector y,
                   NumericMatrix S,
                   double lambda,
                   bool incl_pen = true
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
  
  if(incl_pen){
    // Add penalty term t(beta) %*% S %*% beta
    double pen = 0.0;
    for (int i = 0; i < p; i++) {
      double rowsum = 0.0;
      for(int j = 0; j < p; j ++){
        rowsum +=  S(i,j)  * Betas[j];
      }
      pen += Betas[i] * rowsum;
    }
    sum = sum + pen * 0.5 * lambda;
  }
  
  return(sum);
}