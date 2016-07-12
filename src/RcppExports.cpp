// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// trimFront_cpp
Rcpp::List trimFront_cpp(CharacterVector r_sread, CharacterVector r_qual, CharacterVector r_primer, std::vector<int> prefix_lens);
RcppExport SEXP MotifBinner2_trimFront_cpp(SEXP r_sreadSEXP, SEXP r_qualSEXP, SEXP r_primerSEXP, SEXP prefix_lensSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< CharacterVector >::type r_sread(r_sreadSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_qual(r_qualSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_primer(r_primerSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type prefix_lens(prefix_lensSEXP);
    __result = Rcpp::wrap(trimFront_cpp(r_sread, r_qual, r_primer, prefix_lens));
    return __result;
END_RCPP
}
