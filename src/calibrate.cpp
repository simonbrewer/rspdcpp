#include <Rcpp.h>
using namespace Rcpp;

// Uses interpolation functions taken from
// https://github.com/RfastOfficial/Rfast/issues/12

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
    yout[i] = cpplinterp_(x,y,xout[i],maxgap,false);
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
        while (i-idx>0 && Rcpp::NumericVector::is_na(yout[i-idx]) && x[i]-x[i-idx] < maxgap){
          idx++ ; 
        }
        yout[i]=yout[i-idx];       
      }
    } 
  }
  return yout;
}

// [[Rcpp::export]]
double approx_test(Rcpp::DataFrame calcurve, 
                   int start_date,
                   int end_date){
  
  Rcpp::NumericVector calbp = calcurve["CALBP"];
  Rcpp::NumericVector c14bp = calcurve["C14BP"];
  Rcpp::NumericVector tau1 = calcurve["Error"];
  
  // Interpolation grid
  double xout = 100.0;
  double yout = cpplinterp_(calbp, c14bp, xout, 99999, false);
  return yout;
}

// [[Rcpp::export]]
NumericVector calibrate(Rcpp::NumericVector ages, 
                        Rcpp::NumericVector error, 
                        Rcpp::DataFrame calcurve, 
                        int start_date,
                        int end_date){
  
  // Checks
  if(end_date <= start_date) { 
    Rcpp::stop("end_date should be larger than start_date");
  }
  if(ages.length() != error.length()) { 
    Rcpp::stop("Number of dates does not equal number of errors");
  }
  
  // This is the threshold for inclusion  
  double eps = 1e-5;
  // Number of dates passed
  int n = ages.length();
  
  // Step 1: get yearly conversion series
  Rcpp::NumericVector calbp = calcurve["CALBP"];
  Rcpp::NumericVector c14bp = calcurve["C14BP"];
  Rcpp::NumericVector tau1 = calcurve["Error"];
  
  // Step 2: Make interpolation grid
  IntegerVector tmprange = Rcpp::seq(start_date, end_date);
  NumericVector agegrid = as<NumericVector>(tmprange);
  
  // Step 3: Interpolate Cal and C14 references
  NumericVector mu = cpplinterp(calbp, c14bp, agegrid, 99999, false);
  NumericVector tau2 = cpplinterp(calbp, tau1, agegrid, 99999, false);
  tau2 = pow(tau2, 2);
  
  // Output stores
  NumericVector dens (mu.length());
  NumericVector tau (mu.length());
  
  // Step 4: Convert dates
  for(unsigned int i = 0; i < n; i++) {
    tau = pow(error[i], 2.0) + tau2;
    tau = pow(tau, 0.5);
    dens = cal_dnorm(ages[i], mu, tau, eps);
  }
  return dens;
}

