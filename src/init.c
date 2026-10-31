
// clang-format sorts includes unless SortIncludes: Never. However, the ordering
// does matter here. So, we need to disable clang-format for safety.

// clang-format off
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
// clang-format on

#include "rust/api.h"

static uintptr_t TAGGED_POINTER_MASK = (uintptr_t)1;

SEXP handle_result(SEXP res_) {
    uintptr_t res = (uintptr_t)res_;

    // An error is indicated by tag.
    if ((res & TAGGED_POINTER_MASK) == 1) {
        // Remove tag
        SEXP res_aligned = (SEXP)(res & ~TAGGED_POINTER_MASK);

        // Currently, there are two types of error cases:
        //
        //   1. Error from Rust code
        //   2. Error from R's C API, which is caught by R_UnwindProtect()
        //
        if (TYPEOF(res_aligned) == CHARSXP) {
            // In case 1, the result is an error message that can be passed to
            // Rf_errorcall() directly.
            Rf_errorcall(R_NilValue, "%s", CHAR(res_aligned));
        } else {
            // In case 2, the result is the token to restart the
            // cleanup process on R's side.
            R_ContinueUnwind(res_aligned);
        }
    }

    return (SEXP)res;
}

SEXP savvy_psp_rs_cut_kernel__impl(SEXP c_arg__y, SEXP c_arg__mu, SEXP c_arg__z, SEXP c_arg__kappa_free, SEXP c_arg__n_cat, SEXP c_arg__kmean, SEXP c_arg__ksd, SEXP c_arg__step, SEXP c_arg__reps, SEXP c_arg__ac, SEXP c_arg__mv, SEXP c_arg__seed) {
    SEXP res = savvy_psp_rs_cut_kernel__ffi(c_arg__y, c_arg__mu, c_arg__z, c_arg__kappa_free, c_arg__n_cat, c_arg__kmean, c_arg__ksd, c_arg__step, c_arg__reps, c_arg__ac, c_arg__mv, c_arg__seed);
    return handle_result(res);
}

SEXP savvy_psp_rs_dummy_job__impl(SEXP c_arg__tasks, SEXP c_arg__threads, SEXP c_arg__iters, SEXP c_arg__seed, SEXP c_arg__progress) {
    SEXP res = savvy_psp_rs_dummy_job__ffi(c_arg__tasks, c_arg__threads, c_arg__iters, c_arg__seed, c_arg__progress);
    return handle_result(res);
}

SEXP savvy_psp_rs_filter__impl(SEXP c_arg__w, SEXP c_arg__m, SEXP c_arg__has_m, SEXP c_arg__rho, SEXP c_arg__lambda, SEXP c_arg__t) {
    SEXP res = savvy_psp_rs_filter__ffi(c_arg__w, c_arg__m, c_arg__has_m, c_arg__rho, c_arg__lambda, c_arg__t);
    return handle_result(res);
}

SEXP savvy_psp_rs_ghk__impl(SEXP c_arg__model, SEXP c_arg__draws, SEXP c_arg__y, SEXP c_arg__reps, SEXP c_arg__seed) {
    SEXP res = savvy_psp_rs_ghk__ffi(c_arg__model, c_arg__draws, c_arg__y, c_arg__reps, c_arg__seed);
    return handle_result(res);
}

SEXP savvy_psp_rs_impacts__impl(SEXP c_arg__model, SEXP c_arg__draws, SEXP c_arg__vars_beta, SEXP c_arg__vars_theta, SEXP c_arg__by_period, SEXP c_arg__lesage) {
    SEXP res = savvy_psp_rs_impacts__ffi(c_arg__model, c_arg__draws, c_arg__vars_beta, c_arg__vars_theta, c_arg__by_period, c_arg__lesage);
    return handle_result(res);
}

SEXP savvy_psp_rs_lndet_eigen__impl(SEXP c_arg__re, SEXP c_arg__im, SEXP c_arg__grid) {
    SEXP res = savvy_psp_rs_lndet_eigen__ffi(c_arg__re, c_arg__im, c_arg__grid);
    return handle_result(res);
}

SEXP savvy_psp_rs_lndet_sparse__impl(SEXP c_arg__w, SEXP c_arg__points) {
    SEXP res = savvy_psp_rs_lndet_sparse__ffi(c_arg__w, c_arg__points);
    return handle_result(res);
}

SEXP savvy_psp_rs_next_stream__impl(SEXP c_arg__seed, SEXP c_arg__substream) {
    SEXP res = savvy_psp_rs_next_stream__ffi(c_arg__seed, c_arg__substream);
    return handle_result(res);
}

SEXP savvy_psp_rs_norm_fun__impl(SEXP c_arg__fun, SEXP c_arg__x, SEXP c_arg__lower_tail, SEXP c_arg__log) {
    SEXP res = savvy_psp_rs_norm_fun__ffi(c_arg__fun, c_arg__x, c_arg__lower_tail, c_arg__log);
    return handle_result(res);
}

SEXP savvy_psp_rs_ppc__impl(SEXP c_arg__model, SEXP c_arg__draws, SEXP c_arg__seeds) {
    SEXP res = savvy_psp_rs_ppc__ffi(c_arg__model, c_arg__draws, c_arg__seeds);
    return handle_result(res);
}

SEXP savvy_psp_rs_predict__impl(SEXP c_arg__model, SEXP c_arg__draws) {
    SEXP res = savvy_psp_rs_predict__ffi(c_arg__model, c_arg__draws);
    return handle_result(res);
}

SEXP savvy_psp_rs_probs__impl(SEXP c_arg__model, SEXP c_arg__draws) {
    SEXP res = savvy_psp_rs_probs__ffi(c_arg__model, c_arg__draws);
    return handle_result(res);
}

SEXP savvy_psp_rs_reduced_form__impl(SEXP c_arg__model, SEXP c_arg__draws, SEXP c_arg__need_impacts) {
    SEXP res = savvy_psp_rs_reduced_form__ffi(c_arg__model, c_arg__draws, c_arg__need_impacts);
    return handle_result(res);
}

SEXP savvy_psp_rs_response__impl(SEXP c_arg__model, SEXP c_arg__draws, SEXP c_arg__seeds, SEXP c_arg__spec) {
    SEXP res = savvy_psp_rs_response__ffi(c_arg__model, c_arg__draws, c_arg__seeds, c_arg__spec);
    return handle_result(res);
}

SEXP savvy_psp_rs_rng__impl(SEXP c_arg__kind, SEXP c_arg__n, SEXP c_arg__par, SEXP c_arg__seed) {
    SEXP res = savvy_psp_rs_rng__ffi(c_arg__kind, c_arg__n, c_arg__par, c_arg__seed);
    return handle_result(res);
}

SEXP savvy_psp_rs_rtnorm__impl(SEXP c_arg__n, SEXP c_arg__mean, SEXP c_arg__sd, SEXP c_arg__lo, SEXP c_arg__hi, SEXP c_arg__seed) {
    SEXP res = savvy_psp_rs_rtnorm__ffi(c_arg__n, c_arg__mean, c_arg__sd, c_arg__lo, c_arg__hi, c_arg__seed);
    return handle_result(res);
}

SEXP savvy_psp_rs_sample__impl(SEXP c_arg__model, SEXP c_arg__spec, SEXP c_arg__seeds, SEXP c_arg__threads, SEXP c_arg__progress) {
    SEXP res = savvy_psp_rs_sample__ffi(c_arg__model, c_arg__spec, c_arg__seeds, c_arg__threads, c_arg__progress);
    return handle_result(res);
}

SEXP savvy_psp_rs_weights_apply__impl(SEXP c_arg__w, SEXP c_arg__v) {
    SEXP res = savvy_psp_rs_weights_apply__ffi(c_arg__w, c_arg__v);
    return handle_result(res);
}


static const R_CallMethodDef CallEntries[] = {
    {"savvy_psp_rs_cut_kernel__impl", (DL_FUNC) &savvy_psp_rs_cut_kernel__impl, 12},
    {"savvy_psp_rs_dummy_job__impl", (DL_FUNC) &savvy_psp_rs_dummy_job__impl, 5},
    {"savvy_psp_rs_filter__impl", (DL_FUNC) &savvy_psp_rs_filter__impl, 6},
    {"savvy_psp_rs_ghk__impl", (DL_FUNC) &savvy_psp_rs_ghk__impl, 5},
    {"savvy_psp_rs_impacts__impl", (DL_FUNC) &savvy_psp_rs_impacts__impl, 6},
    {"savvy_psp_rs_lndet_eigen__impl", (DL_FUNC) &savvy_psp_rs_lndet_eigen__impl, 3},
    {"savvy_psp_rs_lndet_sparse__impl", (DL_FUNC) &savvy_psp_rs_lndet_sparse__impl, 2},
    {"savvy_psp_rs_next_stream__impl", (DL_FUNC) &savvy_psp_rs_next_stream__impl, 2},
    {"savvy_psp_rs_norm_fun__impl", (DL_FUNC) &savvy_psp_rs_norm_fun__impl, 4},
    {"savvy_psp_rs_ppc__impl", (DL_FUNC) &savvy_psp_rs_ppc__impl, 3},
    {"savvy_psp_rs_predict__impl", (DL_FUNC) &savvy_psp_rs_predict__impl, 2},
    {"savvy_psp_rs_probs__impl", (DL_FUNC) &savvy_psp_rs_probs__impl, 2},
    {"savvy_psp_rs_reduced_form__impl", (DL_FUNC) &savvy_psp_rs_reduced_form__impl, 3},
    {"savvy_psp_rs_response__impl", (DL_FUNC) &savvy_psp_rs_response__impl, 4},
    {"savvy_psp_rs_rng__impl", (DL_FUNC) &savvy_psp_rs_rng__impl, 4},
    {"savvy_psp_rs_rtnorm__impl", (DL_FUNC) &savvy_psp_rs_rtnorm__impl, 6},
    {"savvy_psp_rs_sample__impl", (DL_FUNC) &savvy_psp_rs_sample__impl, 5},
    {"savvy_psp_rs_weights_apply__impl", (DL_FUNC) &savvy_psp_rs_weights_apply__impl, 2},
    {NULL, NULL, 0}
};

void R_init_panelspatialprobit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    // Functions for initialization, if any.

}
