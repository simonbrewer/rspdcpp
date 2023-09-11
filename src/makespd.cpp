#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List  make_spd(Rcpp::List input_val, 
                             int start_date,
                             int end_date){
  
  if(end_date <= start_date) { 
    Rcpp::stop("end_date should be larger than start_date");
  }
  // Create vector to store sum of probability
  Rcpp::NumericVector prsum ((end_date + 1) - start_date);
  Rcpp::IntegerVector ages = Rcpp::seq(start_date, end_date);
  
  // Get size of lists (no of dates)
  unsigned int n = input_val.length();
  unsigned int elems = 0;
  
  // loop over dates
  for(unsigned int i = 0; i < n; i++) {
    Rcpp::DataFrame df(input_val[i]);
    Rcpp::IntegerVector calbp = df["calBP"];
    Rcpp::NumericVector prdens = df["PrDens"];
    //Rcpp::DataFrame dates = input_val[i];
    elems = calbp.length();
    for(unsigned int j = 0; j < elems; j++) {
      if (calbp(j) >= start_date && calbp(j) <= end_date) {
        prsum(calbp(j)) += prdens(j);
      }
    }
    
  }
  prsum = prsum / sum(prsum);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("calbp") = ages,
                                      Rcpp::Named("d") = prsum);
  return out;
}

