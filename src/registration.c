/* Native routine registration while the package holds both the C++ (Rcpp)
 * engine and the Rust code (plan/06, phases R0-R4).
 *
 * savvy's generated init.c is compiled as part of this file (src/Makevars
 * lists OBJECTS without init.o), with its R_init renamed. Its CallEntries
 * table is copied at load time, so new savvy functions need no edit here.
 * Rcpp::compileAttributes() sees the R_init below and leaves registration to
 * it; the frozen Rcpp routines are listed by hand.
 * Delete this file in R5, when init.c becomes the only registration. */

#define R_init_panelspatialprobit psp_savvy_init
#include "init.c"
#undef R_init_panelspatialprobit

#include <string.h>
#include <R_ext/Rdynload.h>

/* Rcpp: src/RcppExports.cpp */
SEXP _panelspatialprobit_psp_ghk_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _panelspatialprobit_psp_sampler_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                         SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                         SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                         SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                         SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                         SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _panelspatialprobit_psp_rtnorm_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _panelspatialprobit_psp_cut_kernel_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                            SEXP, SEXP, SEXP);

static const R_CallMethodDef RcppEntries[] = {
    {"_panelspatialprobit_psp_ghk_cpp", (DL_FUNC) &_panelspatialprobit_psp_ghk_cpp, 6},
    {"_panelspatialprobit_psp_sampler_cpp", (DL_FUNC) &_panelspatialprobit_psp_sampler_cpp, 50},
    {"_panelspatialprobit_psp_rtnorm_cpp", (DL_FUNC) &_panelspatialprobit_psp_rtnorm_cpp, 5},
    {"_panelspatialprobit_psp_cut_kernel_cpp", (DL_FUNC) &_panelspatialprobit_psp_cut_kernel_cpp, 11},
    {NULL, NULL, 0}
};

#define N_RCPP (sizeof(RcppEntries) / sizeof(RcppEntries[0]) - 1)
#define N_SAVVY (sizeof(CallEntries) / sizeof(CallEntries[0]) - 1)
static R_CallMethodDef AllEntries[N_RCPP + N_SAVVY + 1];

void R_init_panelspatialprobit(DllInfo *dll) {
    /* savvy's own init first, for any initialization calls it makes; the
     * second registration replaces its routine table with the merged one. */
    psp_savvy_init(dll);
    memcpy(AllEntries, RcppEntries, N_RCPP * sizeof(R_CallMethodDef));
    memcpy(AllEntries + N_RCPP, CallEntries, (N_SAVVY + 1) * sizeof(R_CallMethodDef));
    R_registerRoutines(dll, NULL, AllEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
