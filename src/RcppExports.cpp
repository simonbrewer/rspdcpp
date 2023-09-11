// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cal_dnorm
Rcpp::NumericVector cal_dnorm(double x, Rcpp::NumericVector means, Rcpp::NumericVector sds, double eps);
RcppExport SEXP _rspdcpp_cal_dnorm(SEXP xSEXP, SEXP meansSEXP, SEXP sdsSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type means(meansSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sds(sdsSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_dnorm(x, means, sds, eps));
    return rcpp_result_gen;
END_RCPP
}
// cpplinterp_
double cpplinterp_(Rcpp::NumericVector x, Rcpp::NumericVector y, double xout, const int maxgap, const bool locf);
RcppExport SEXP _rspdcpp_cpplinterp_(SEXP xSEXP, SEXP ySEXP, SEXP xoutSEXP, SEXP maxgapSEXP, SEXP locfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< const int >::type maxgap(maxgapSEXP);
    Rcpp::traits::input_parameter< const bool >::type locf(locfSEXP);
    rcpp_result_gen = Rcpp::wrap(cpplinterp_(x, y, xout, maxgap, locf));
    return rcpp_result_gen;
END_RCPP
}
// cpplinterp
Rcpp::NumericVector cpplinterp(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector xout, const int maxgap, const bool locf);
RcppExport SEXP _rspdcpp_cpplinterp(SEXP xSEXP, SEXP ySEXP, SEXP xoutSEXP, SEXP maxgapSEXP, SEXP locfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< const int >::type maxgap(maxgapSEXP);
    Rcpp::traits::input_parameter< const bool >::type locf(locfSEXP);
    rcpp_result_gen = Rcpp::wrap(cpplinterp(x, y, xout, maxgap, locf));
    return rcpp_result_gen;
END_RCPP
}
// approx_test
double approx_test(Rcpp::DataFrame calcurve, int start_date, int end_date);
RcppExport SEXP _rspdcpp_approx_test(SEXP calcurveSEXP, SEXP start_dateSEXP, SEXP end_dateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type calcurve(calcurveSEXP);
    Rcpp::traits::input_parameter< int >::type start_date(start_dateSEXP);
    Rcpp::traits::input_parameter< int >::type end_date(end_dateSEXP);
    rcpp_result_gen = Rcpp::wrap(approx_test(calcurve, start_date, end_date));
    return rcpp_result_gen;
END_RCPP
}
// calibrate
List calibrate(Rcpp::NumericVector ages, Rcpp::NumericVector error, Rcpp::DataFrame calcurve, int start_date, int end_date, const bool normalize);
RcppExport SEXP _rspdcpp_calibrate(SEXP agesSEXP, SEXP errorSEXP, SEXP calcurveSEXP, SEXP start_dateSEXP, SEXP end_dateSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type error(errorSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type calcurve(calcurveSEXP);
    Rcpp::traits::input_parameter< int >::type start_date(start_dateSEXP);
    Rcpp::traits::input_parameter< int >::type end_date(end_dateSEXP);
    Rcpp::traits::input_parameter< const bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(calibrate(ages, error, calcurve, start_date, end_date, normalize));
    return rcpp_result_gen;
END_RCPP
}
// uncalibrate
List uncalibrate(Rcpp::NumericVector agegrid, Rcpp::NumericVector prdens, Rcpp::DataFrame calcurve, int start_date, int end_date, const bool normalize);
RcppExport SEXP _rspdcpp_uncalibrate(SEXP agegridSEXP, SEXP prdensSEXP, SEXP calcurveSEXP, SEXP start_dateSEXP, SEXP end_dateSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type agegrid(agegridSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prdens(prdensSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type calcurve(calcurveSEXP);
    Rcpp::traits::input_parameter< int >::type start_date(start_dateSEXP);
    Rcpp::traits::input_parameter< int >::type end_date(end_dateSEXP);
    Rcpp::traits::input_parameter< const bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(uncalibrate(agegrid, prdens, calcurve, start_date, end_date, normalize));
    return rcpp_result_gen;
END_RCPP
}
// make_spd
Rcpp::List make_spd(Rcpp::List input_val, int start_date, int end_date);
RcppExport SEXP _rspdcpp_make_spd(SEXP input_valSEXP, SEXP start_dateSEXP, SEXP end_dateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type input_val(input_valSEXP);
    Rcpp::traits::input_parameter< int >::type start_date(start_dateSEXP);
    Rcpp::traits::input_parameter< int >::type end_date(end_dateSEXP);
    rcpp_result_gen = Rcpp::wrap(make_spd(input_val, start_date, end_date));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rspdcpp_cal_dnorm", (DL_FUNC) &_rspdcpp_cal_dnorm, 4},
    {"_rspdcpp_cpplinterp_", (DL_FUNC) &_rspdcpp_cpplinterp_, 5},
    {"_rspdcpp_cpplinterp", (DL_FUNC) &_rspdcpp_cpplinterp, 5},
    {"_rspdcpp_approx_test", (DL_FUNC) &_rspdcpp_approx_test, 3},
    {"_rspdcpp_calibrate", (DL_FUNC) &_rspdcpp_calibrate, 6},
    {"_rspdcpp_uncalibrate", (DL_FUNC) &_rspdcpp_uncalibrate, 6},
    {"_rspdcpp_make_spd", (DL_FUNC) &_rspdcpp_make_spd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_rspdcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
