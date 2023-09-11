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
                   double xout, const int maxgap=99999, 
                   const bool locf=false){
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
List calibrate(Rcpp::NumericVector ages, 
               Rcpp::NumericVector error, 
               Rcpp::DataFrame calcurve, 
               int start_date,
               int end_date,
               const bool normalize=false){
  
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
  // Output list
  Rcpp::List output (n);
  
  
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
    if (normalize == true) {
      dens = dens / sum(dens);
    }
    
    Rcpp::DataFrame tmp =
      Rcpp::DataFrame::create(Rcpp::Named("calbp") = agegrid[dens > eps],
                              Rcpp::Named("prdens") = dens[dens > eps]);
    output[i] = tmp;
  }
  return output;
}

// [[Rcpp::export]]
List uncalibrate(Rcpp::NumericVector agegrid, 
                 Rcpp::NumericVector prdens, 
                 Rcpp::DataFrame calcurve, 
                 int start_date,
                 int end_date,
                 const bool normalize=false){
  
  // Checks
  if(end_date <= start_date) { 
    Rcpp::stop("end_date should be larger than start_date");
  }
  if(agegrid.length() != prdens.length()) { 
    Rcpp::stop("Number of dates does not equal number of densities");
  }
  
  // This is the threshold for inclusion  
  double eps = 0;
  // Number of dates passed
  int n = agegrid.length();
  
  // Step 1: get yearly conversion series
  Rcpp::NumericVector calbp = calcurve["CALBP"];
  Rcpp::NumericVector c14bp = calcurve["C14BP"];
  Rcpp::NumericVector error = calcurve["Error"];
  
  // Step 2: Make interpolation grid
  
  // Step 3: Interpolate Cal and C14 references
  NumericVector c14grid = cpplinterp(calbp, c14bp, agegrid, 99999, false);
  NumericVector c14err = cpplinterp(calbp, error, agegrid, 99999, false);
  
  NumericVector rCRA (c14grid.length());
  for(unsigned int i = 0; i < n; i++) {
    rCRA[i] = round(R::rnorm(c14grid[i], c14err[i]));
  }
  
  Rcpp::NumericVector h = prdens/sum(prdens);
  // mu = mycras$ccCRA; THIS IS C14GRID
  //  s = mycras$ccError: THIS IS C14ERR
  Rcpp::IntegerVector k = Rcpp::seq(start_date, end_date);
  Rcpp::NumericVector base (k.length());
  Rcpp::NumericVector raw (k.length());
  Rcpp::NumericVector dens (k.length());
  Rcpp::NumericVector tmp (c14grid.length());
  
  for(unsigned int i = 0; i < k.length(); i++) {
    tmp = cal_dnorm(k[i], c14grid, c14err, eps);
    base[i] = sum(tmp);
    raw[i] = sum(tmp * h);
  }
  
  // Normalizing (i think)
  raw = raw/sum(raw);
  raw[raw < 1e-5] = 0;
  dens[base > 0] = raw[base > 0] / base[base > 0];
  dens = dens / sum(dens);
  // for(unsigned int i = 0; i < k.length(); i++) {
  //   
  // }
  // Output
  Rcpp::List output =
    Rcpp::List::create(Rcpp::Named("cra") = k[dens > 0],
                       Rcpp::Named("raw") = raw[dens > 0],
                                               Rcpp::Named("base") = base[dens > 0],
                                                                         Rcpp::Named("prdens") = dens[dens > 0]);
  
  return output;
}

