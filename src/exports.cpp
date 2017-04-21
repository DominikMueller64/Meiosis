#include "../inst/include/Meiosis.h"
#include <Rcpp.h>

using namespace Rcpp;

// seed_rng
int seed_rng(const Rcpp::Nullable<int>& seed);
RcppExport SEXP Meiosis_seed_rng(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::Nullable<int>& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(seed_rng(seed));
    return rcpp_result_gen;
END_RCPP
}
// calc_Lstar
double calc_Lstar(const double L, const int m, const double p, const Rcpp::Nullable<double> epsilon);
RcppExport SEXP Meiosis_calc_Lstar(SEXP LSEXP, SEXP mSEXP, SEXP pSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<double> >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_Lstar(L, m, p, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// crossover
Rcpp::NumericVector crossover(const double L, const int m, const double p, const bool obligate_chiasma, const double Lstar);
RcppExport SEXP Meiosis_crossover(SEXP LSEXP, SEXP mSEXP, SEXP pSEXP, SEXP obligate_chiasmaSEXP, SEXP LstarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const bool >::type obligate_chiasma(obligate_chiasmaSEXP);
    Rcpp::traits::input_parameter< const double >::type Lstar(LstarSEXP);
    rcpp_result_gen = Rcpp::wrap(crossover(L, m, p, obligate_chiasma, Lstar));
    return rcpp_result_gen;
END_RCPP
}
// meiosis_geno
Rcpp::List meiosis_geno(const Rcpp::List& individual, const Rcpp::List& positions, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_meiosis_geno(SEXP individualSEXP, SEXP positionsSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(meiosis_geno(individual, positions, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// meiosis_xo
Rcpp::List meiosis_xo(const Rcpp::List& individual, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_meiosis_xo(SEXP individualSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(meiosis_xo(individual, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// cross_geno
Rcpp::List cross_geno(const Rcpp::List& father, const Rcpp::List& mother, const Rcpp::List& positions, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_cross_geno(SEXP fatherSEXP, SEXP motherSEXP, SEXP positionsSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type father(fatherSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_geno(father, mother, positions, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// cross_xo
Rcpp::List cross_xo(const Rcpp::List& father, const Rcpp::List& mother, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_cross_xo(SEXP fatherSEXP, SEXP motherSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type father(fatherSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_xo(father, mother, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// self_geno
Rcpp::List self_geno(const Rcpp::List& individual, const Rcpp::List& positions, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_self_geno(SEXP individualSEXP, SEXP positionsSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(self_geno(individual, positions, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// self_xo
Rcpp::List self_xo(const Rcpp::List& individual, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_self_xo(SEXP individualSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(self_xo(individual, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// dh_geno
Rcpp::List dh_geno(const Rcpp::List& individual, const Rcpp::List& positions, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_dh_geno(SEXP individualSEXP, SEXP positionsSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(dh_geno(individual, positions, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// dh_xo
Rcpp::List dh_xo(const Rcpp::List& individual, const Rcpp::List& xoparam);
RcppExport SEXP Meiosis_dh_xo(SEXP individualSEXP, SEXP xoparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xoparam(xoparamSEXP);
    rcpp_result_gen = Rcpp::wrap(dh_xo(individual, xoparam));
    return rcpp_result_gen;
END_RCPP
}
// realized_heter
double realized_heter(const Rcpp::List& individual);
RcppExport SEXP Meiosis_realized_heter(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(realized_heter(individual));
    return rcpp_result_gen;
END_RCPP
}
// realized_coancestry
double realized_coancestry(const Rcpp::List& individual_1, const Rcpp::Nullable<Rcpp::List>& individual_2);
RcppExport SEXP Meiosis_realized_coancestry(SEXP individual_1SEXP, SEXP individual_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual_1(individual_1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::List>& >::type individual_2(individual_2SEXP);
    rcpp_result_gen = Rcpp::wrap(realized_coancestry(individual_1, individual_2));
    return rcpp_result_gen;
END_RCPP
}
