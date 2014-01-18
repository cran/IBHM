#include <Rcpp.h>

using namespace Rcpp;

// weightedMean
double weightedMean(const NumericVector& x, const NumericVector& w);
RcppExport SEXP IBHM_weightedMean(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    NumericVector x = Rcpp::as<NumericVector >(xSEXP);
    NumericVector w = Rcpp::as<NumericVector >(wSEXP);
    double __result = weightedMean(x, w);
    return Rcpp::wrap(__result);
END_RCPP
}
// weightedVar
double weightedVar(const NumericVector& x, const NumericVector& w);
RcppExport SEXP IBHM_weightedVar(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    NumericVector x = Rcpp::as<NumericVector >(xSEXP);
    NumericVector w = Rcpp::as<NumericVector >(wSEXP);
    double __result = weightedVar(x, w);
    return Rcpp::wrap(__result);
END_RCPP
}
// weightedCov
double weightedCov(const NumericVector& x, const NumericVector& z, const NumericVector& w);
RcppExport SEXP IBHM_weightedCov(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP) {
BEGIN_RCPP
    NumericVector x = Rcpp::as<NumericVector >(xSEXP);
    NumericVector z = Rcpp::as<NumericVector >(zSEXP);
    NumericVector w = Rcpp::as<NumericVector >(wSEXP);
    double __result = weightedCov(x, z, w);
    return Rcpp::wrap(__result);
END_RCPP
}
// weightedR
double weightedR(const NumericVector& y1, const NumericVector& y2, const NumericVector& w);
RcppExport SEXP IBHM_weightedR(SEXP y1SEXP, SEXP y2SEXP, SEXP wSEXP) {
BEGIN_RCPP
    NumericVector y1 = Rcpp::as<NumericVector >(y1SEXP);
    NumericVector y2 = Rcpp::as<NumericVector >(y2SEXP);
    NumericVector w = Rcpp::as<NumericVector >(wSEXP);
    double __result = weightedR(y1, y2, w);
    return Rcpp::wrap(__result);
END_RCPP
}
// tiedRanks
NumericVector tiedRanks(const NumericVector& x);
RcppExport SEXP IBHM_tiedRanks(SEXP xSEXP) {
BEGIN_RCPP
    NumericVector x = Rcpp::as<NumericVector >(xSEXP);
    NumericVector __result = tiedRanks(x);
    return Rcpp::wrap(__result);
END_RCPP
}
// weightedRho
double weightedRho(const NumericVector& y1, const NumericVector& y2, const NumericVector& w);
RcppExport SEXP IBHM_weightedRho(SEXP y1SEXP, SEXP y2SEXP, SEXP wSEXP) {
BEGIN_RCPP
    NumericVector y1 = Rcpp::as<NumericVector >(y1SEXP);
    NumericVector y2 = Rcpp::as<NumericVector >(y2SEXP);
    NumericVector w = Rcpp::as<NumericVector >(wSEXP);
    double __result = weightedRho(y1, y2, w);
    return Rcpp::wrap(__result);
END_RCPP
}
