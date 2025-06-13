// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// blockwise_crossprod
arma::mat blockwise_crossprod(const arma::mat& X, int n_threads, int block_size);
RcppExport SEXP _SuSiE4X_blockwise_crossprod(SEXP XSEXP, SEXP n_threadsSEXP, SEXP block_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type block_size(block_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(blockwise_crossprod(X, n_threads, block_size));
    return rcpp_result_gen;
END_RCPP
}
// blockwise_crossprod2
arma::mat blockwise_crossprod2(const arma::mat& X, const arma::mat& Z, int n_threads, int block_size);
RcppExport SEXP _SuSiE4X_blockwise_crossprod2(SEXP XSEXP, SEXP ZSEXP, SEXP n_threadsSEXP, SEXP block_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type block_size(block_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(blockwise_crossprod2(X, Z, n_threads, block_size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SuSiE4X_blockwise_crossprod", (DL_FUNC) &_SuSiE4X_blockwise_crossprod, 3},
    {"_SuSiE4X_blockwise_crossprod2", (DL_FUNC) &_SuSiE4X_blockwise_crossprod2, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SuSiE4X(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
