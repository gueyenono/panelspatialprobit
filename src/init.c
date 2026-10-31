
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

SEXP savvy_psp_rs_dummy_job__impl(SEXP c_arg__tasks, SEXP c_arg__threads, SEXP c_arg__iters, SEXP c_arg__seed, SEXP c_arg__progress) {
    SEXP res = savvy_psp_rs_dummy_job__ffi(c_arg__tasks, c_arg__threads, c_arg__iters, c_arg__seed, c_arg__progress);
    return handle_result(res);
}

SEXP savvy_psp_rs_filter__impl(SEXP c_arg__w, SEXP c_arg__m, SEXP c_arg__has_m, SEXP c_arg__rho, SEXP c_arg__lambda, SEXP c_arg__t) {
    SEXP res = savvy_psp_rs_filter__ffi(c_arg__w, c_arg__m, c_arg__has_m, c_arg__rho, c_arg__lambda, c_arg__t);
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

SEXP savvy_psp_rs_rng__impl(SEXP c_arg__kind, SEXP c_arg__n, SEXP c_arg__par, SEXP c_arg__seed) {
    SEXP res = savvy_psp_rs_rng__ffi(c_arg__kind, c_arg__n, c_arg__par, c_arg__seed);
    return handle_result(res);
}

SEXP savvy_psp_rs_rtnorm__impl(SEXP c_arg__n, SEXP c_arg__mean, SEXP c_arg__sd, SEXP c_arg__lo, SEXP c_arg__hi, SEXP c_arg__seed) {
    SEXP res = savvy_psp_rs_rtnorm__ffi(c_arg__n, c_arg__mean, c_arg__sd, c_arg__lo, c_arg__hi, c_arg__seed);
    return handle_result(res);
}

SEXP savvy_psp_rs_weights_apply__impl(SEXP c_arg__w, SEXP c_arg__v) {
    SEXP res = savvy_psp_rs_weights_apply__ffi(c_arg__w, c_arg__v);
    return handle_result(res);
}


static const R_CallMethodDef CallEntries[] = {
    {"savvy_psp_rs_dummy_job__impl", (DL_FUNC) &savvy_psp_rs_dummy_job__impl, 5},
    {"savvy_psp_rs_filter__impl", (DL_FUNC) &savvy_psp_rs_filter__impl, 6},
    {"savvy_psp_rs_next_stream__impl", (DL_FUNC) &savvy_psp_rs_next_stream__impl, 2},
    {"savvy_psp_rs_norm_fun__impl", (DL_FUNC) &savvy_psp_rs_norm_fun__impl, 4},
    {"savvy_psp_rs_rng__impl", (DL_FUNC) &savvy_psp_rs_rng__impl, 4},
    {"savvy_psp_rs_rtnorm__impl", (DL_FUNC) &savvy_psp_rs_rtnorm__impl, 6},
    {"savvy_psp_rs_weights_apply__impl", (DL_FUNC) &savvy_psp_rs_weights_apply__impl, 2},
    {NULL, NULL, 0}
};

void R_init_panelspatialprobit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    // Functions for initialization, if any.

}
