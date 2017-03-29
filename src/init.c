#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP Meiosis_bcgv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Meiosis_calc_Lstar(SEXP, SEXP, SEXP, SEXP);
extern SEXP Meiosis_cross_geno(SEXP, SEXP, SEXP, SEXP);
extern SEXP Meiosis_cross_xo(SEXP, SEXP, SEXP);
extern SEXP Meiosis_crossover(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Meiosis_dh_geno(SEXP, SEXP, SEXP);
extern SEXP Meiosis_dh_xo(SEXP, SEXP);
extern SEXP Meiosis_meiosis_geno(SEXP, SEXP, SEXP);
extern SEXP Meiosis_meiosis_xo(SEXP, SEXP);
extern SEXP Meiosis_realized_coancestry(SEXP, SEXP);
extern SEXP Meiosis_seed_rng(SEXP);
extern SEXP Meiosis_to_matrix(SEXP);
extern SEXP _rcpp_module_boot_Module(void);

static const R_CallMethodDef CallEntries[] =
{
    {"Meiosis_bcgv",                (DL_FUNC) &Meiosis_bcgv,                9},
    {"Meiosis_calc_Lstar",          (DL_FUNC) &Meiosis_calc_Lstar,          4},
    {"Meiosis_cross_geno",          (DL_FUNC) &Meiosis_cross_geno,          4},
    {"Meiosis_cross_xo",            (DL_FUNC) &Meiosis_cross_xo,            3},
    {"Meiosis_crossover",           (DL_FUNC) &Meiosis_crossover,           5},
    {"Meiosis_dh_geno",             (DL_FUNC) &Meiosis_dh_geno,             3},
    {"Meiosis_dh_xo",               (DL_FUNC) &Meiosis_dh_xo,               2},
    {"Meiosis_meiosis_geno",        (DL_FUNC) &Meiosis_meiosis_geno,        3},
    {"Meiosis_meiosis_xo",          (DL_FUNC) &Meiosis_meiosis_xo,          2},
    {"Meiosis_realized_coancestry", (DL_FUNC) &Meiosis_realized_coancestry, 2},
    {"Meiosis_seed_rng",            (DL_FUNC) &Meiosis_seed_rng,            1},
    {"Meiosis_to_matrix",           (DL_FUNC) &Meiosis_to_matrix,           1},
    {"_rcpp_module_boot_Module",    (DL_FUNC) &_rcpp_module_boot_Module,    0},
    {NULL,                          NULL,                                   0}
};

void R_init_Meiosis(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
