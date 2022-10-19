// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ComputeAlphaBeta
List ComputeAlphaBeta(arma::vec y, arma::vec x, arma::mat WW, arma::vec weightVec, arma::mat Z, arma::vec G);
RcppExport SEXP _bartik_weight_ComputeAlphaBeta(SEXP ySEXP, SEXP xSEXP, SEXP WWSEXP, SEXP weightVecSEXP, SEXP ZSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type WW(WWSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weightVec(weightVecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeAlphaBeta(y, x, WW, weightVec, Z, G));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bartik_weight_ComputeAlphaBeta", (DL_FUNC) &_bartik_weight_ComputeAlphaBeta, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bartik_weight(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
