#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix make_mat(Rcpp::List input_list){
  
  unsigned int n = input_list.length();
  
  if(n == 0) { 
    Rcpp::stop("Must supply a list with more than 1 element.");
  }
  
  Rcpp::NumericVector testvals = input_list[0];
  unsigned int elems = testvals.length();
  
  Rcpp::NumericMatrix result_mat = Rcpp::no_init(n, elems);
  
  // fill by row
  for(unsigned int i = 0; i < n; i++) {
    Rcpp::NumericVector row_val = input_list[i];
    
    if(elems != row_val.length()) {
      Rcpp::stop("Length of row does not match matrix requirements"); 
    }
    
    result_mat(i, Rcpp::_) = row_val;
    
  }
  
  return result_mat;
}

// [[Rcpp::export]]
int print_list(Rcpp::List input_list){
  
  unsigned int n = input_list.length();
  
  if(n == 0) { 
    Rcpp::stop("Must supply a list with more than 1 element.");
  }
  
  Rcpp::NumericVector testvals = input_list[0];
  unsigned int elems = testvals.length();
  
  return elems;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix make_mat_order(Rcpp::List input_val, 
                                   Rcpp::List input_pos,
                                   int start_date,
                                   int end_date){
  
  if(end_date <= start_date) { 
    Rcpp::stop("end_date should be larger than start_date");
  }
  
  // Create vector to store sum of probability
  Rcpp::NumericVector prsum (end_date - start_date);
  
  // Get size of lists (no of dates)
  unsigned int n0 = input_val.length();
  unsigned int n1 = input_pos.length();
  
  if(n0 == 0) { 
    Rcpp::stop("Must supply a list with more than 1 element.");
  }
  if(n0 != n1) { 
    Rcpp::stop("Length of values and positions must be identical");
  }
  
  Rcpp::NumericVector testvals = input_val[0];
  unsigned int elems = testvals.length();
  
  Rcpp::NumericMatrix result_mat = Rcpp::no_init(n0, elems);
  
  // fill by row
  for(unsigned int i = 0; i < n0; i++) {
    Rcpp::NumericVector row_val = input_val[i];
    Rcpp::NumericVector row_pos = input_pos[i];
    
    if(elems != row_val.length()) {
      Rcpp::stop("Length of row does not match matrix requirements"); 
    }
    for(unsigned int j = 0; j < elems; j++) {
      result_mat(i, row_pos(j)) = row_val(j);
      prsum(j) += row_val(j);
    }
    
  }
  
  //return result_mat;
  return result_mat;
}

// [[Rcpp::export]]
Rcpp::NumericVector make_spd(Rcpp::List input_val, 
                             Rcpp::List input_pos,
                             int start_date,
                             int end_date){
  
  if(end_date <= start_date) { 
    Rcpp::stop("end_date should be larger than start_date");
  }
  
  // Create vector to store sum of probability
  Rcpp::NumericVector prsum (end_date - start_date);
  
  // Get size of lists (no of dates)
  unsigned int n0 = input_val.length();
  unsigned int n1 = input_pos.length();
  
  if(n0 == 0) { 
    Rcpp::stop("Must supply a list with more than 1 element.");
  }
  if(n0 != n1) { 
    Rcpp::stop("Length of values and positions must be identical");
  }
  
  Rcpp::NumericVector testvals = input_val[0];
  unsigned int elems = testvals.length();
  
  Rcpp::NumericMatrix result_mat = Rcpp::no_init(n0, elems);
  
  // fill by row
  for(unsigned int i = 0; i < n0; i++) {
    Rcpp::NumericVector row_val = input_val[i];
    Rcpp::NumericVector row_pos = input_pos[i];
    
    if(elems != row_val.length()) {
      Rcpp::stop("Length of row does not match matrix requirements"); 
    }
    for(unsigned int j = 0; j < elems; j++) {
      result_mat(i, row_pos(j)) = row_val(j);
      prsum(j) += row_val(j);
    }
    
  }
  prsum = prsum / n0;
  //return result_mat;
  return prsum;
}

