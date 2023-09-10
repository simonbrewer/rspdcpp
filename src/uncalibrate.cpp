#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector cal_dnorm( double x, 
                               Rcpp::NumericVector means, 
                               Rcpp::NumericVector sds,
                               double eps){
  int n = means.size() ;
  Rcpp::NumericVector res(n) ;
  for( int i=0; i<n; i++){
    res[i] = R::dnorm( x, means[i], sds[i], false );
    if (res[i] < eps) res[i] = 0;
  } 
  return res ;
}


