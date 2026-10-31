/* Native routine registration while the package holds both the C++ (Rcpp)
 * engine and the Rust code (plan/06, phases R0-R4).
 *
 * Rcpp::compileAttributes() sees the R_init function below and leaves the
 * registration to it. savvy's generated init.c defines its own R_init, which
 * src/Makevars renames to psp_savvy_init (it is never called). So:
 *   - after Rcpp::compileAttributes(), copy any changed Rcpp entry from the
 *     declarations in RcppExports.cpp;
 *   - after savvy-cli update, add each new routine from the CallEntries of
 *     init.c ("savvy_<fn>__impl"), and mirror any initialization calls that
 *     savvy adds to its R_init.
 * Delete this file in R5, when init.c becomes the only registration. */

#undef R_init_panelspatialprobit

#include <R.h>
#include <Rinternals.h>
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

/* savvy: src/init.c */
SEXP savvy_psp_rs_rng__impl(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_panelspatialprobit_psp_ghk_cpp", (DL_FUNC) &_panelspatialprobit_psp_ghk_cpp, 6},
    {"_panelspatialprobit_psp_sampler_cpp", (DL_FUNC) &_panelspatialprobit_psp_sampler_cpp, 50},
    {"_panelspatialprobit_psp_rtnorm_cpp", (DL_FUNC) &_panelspatialprobit_psp_rtnorm_cpp, 5},
    {"_panelspatialprobit_psp_cut_kernel_cpp", (DL_FUNC) &_panelspatialprobit_psp_cut_kernel_cpp, 11},
    {"savvy_psp_rs_rng__impl", (DL_FUNC) &savvy_psp_rs_rng__impl, 3},
    {NULL, NULL, 0}
};

void R_init_panelspatialprobit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
