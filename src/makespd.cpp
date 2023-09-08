#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cpplinterp_(Rcpp::NumericVector x,
                   Rcpp::NumericVector y, 
                   double xout, const int maxgap=99999, const bool locf=false){
  // Linear approximation of a vector-function
  // x, y	vectors giving the coordinates of the points to be interpolated. 
  // x is assumed to be strictly monotonic.
  // xout	points at which to interpolate.
  int n_x=x.size();
  if( xout < x[0] ){
    if(locf){
      return y[0];
    }else{
      return NA_REAL;
    }
  }
  if( xout >= x[n_x-1] ){
    if(locf){
      return y[n_x-1];
    }else if(xout>x[n_x-1]){
      return NA_REAL;
    }else{
      return y[n_x-1];
    }    
  }
  
  int idx = 0 ;// Find the location of xout in x
  while ( x[idx] <= xout && idx < n_x ){
    idx++ ; 
  } 
  // Find the nonNA location of y around xout within maxgap
  int leftjust=1;
  int rightjust=0;
  while (Rcpp::NumericVector::is_na(y[idx-leftjust]) && idx-leftjust > 0 ){
    leftjust++; 
  }
  while (Rcpp::NumericVector::is_na(y[idx+rightjust]) && idx+rightjust < n_x ){
    rightjust++ ; 
  }
  if(x[idx+rightjust] - x[idx-leftjust] -1<=maxgap){
    double theta = ( xout - x[idx-leftjust] ) / ( x[idx+rightjust] - x[idx-leftjust] ) ; 
    return theta * y[idx+rightjust] + ( 1 - theta ) * y[idx-leftjust] ;  
  }else{
    return NA_REAL;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector cpplinterp(Rcpp::NumericVector x, 
                         Rcpp::NumericVector y, 
                         Rcpp::NumericVector xout, 
                         const int maxgap=99999, const bool locf=false){
  Rcpp::NumericVector yout=clone(xout);
  for(int i=0;i<yout.size();i++){
    yout[i]=cpplinterp_(x,y,xout[i],maxgap,false);
  }
  if(locf){
    for(int i=0;i<yout.size();i++){
      if(Rcpp::NumericVector::is_na(yout[i])){
        int idx=1;
        while (i+idx< yout.size()-1 && 
               Rcpp::NumericVector::is_na(yout[i+idx]) && 
               x[i+idx]-x[i] < maxgap){
          idx++ ; 
        }
        yout[i]=yout[i+idx];       
      }
    }     
    for(int i=yout.size()-1;i>=0;i--){
      if(Rcpp::NumericVector::is_na(yout[i])){
        int idx=1;
        while (i-idx>0 && NumericVector::is_na(yout[i-idx]) && x[i]-x[i-idx] < maxgap){
          idx++ ; 
        }
        yout[i]=yout[i-idx];       
      }
    } 
  }
  return yout;
}

// [[Rcpp::export]]
NumericVector calibrate(){
  return NA_REAL;
}

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

// [[Rcpp::export]]
Rcpp::List test_spd(Rcpp::List input_val, 
                     int start_date,
                     int end_date){
  
  if(end_date <= start_date) { 
    Rcpp::stop("end_date should be larger than start_date");
  }
  // Create vector to store sum of probability
  Rcpp::NumericVector prsum ((end_date + 1) - start_date);
  Rcpp::IntegerVector ages = Rcpp::seq(start_date, end_date);
  
  Rcpp::DataFrame df(input_val[0]);
  Rcpp::IntegerVector calbp = df["calBP"];
  Rcpp::NumericVector prdens = df["PrDens"];
  prdens = prdens * 1.0;
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("calbp") = calbp,
                                      Rcpp::Named("d") = prdens);
  return out;
}
