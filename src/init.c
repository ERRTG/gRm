#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

/*
 * One registration table is the native API boundary. Keep arithmetic kernels
 * in their source-routine files and declare only their R-callable entry points
 * here. Validation-only solvers must not be registered in the package DLL.
 */
extern SEXP gRm_screen_j_exact_chi_gamma_slices(
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);
extern SEXP gRm_screen_j_exact_chi_gamma_trace_slices(
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);
extern SEXP gRm_screen_j_exact_chi_slices(
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);
extern SEXP gRm_screen_j_exact_gamma_slices(SEXP, SEXP, SEXP, SEXP);
extern SEXP gRm_screen_j_conditional_bias_test(
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);
extern SEXP gRm_screen_j_item_pair_conditional_exact(
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);
extern SEXP gRm_screen_j_source_random_draws(SEXP, SEXP);
extern SEXP gRm_item_parameters_extended_capabilities(void);
extern SEXP gRm_item_parameters_top_ice_field(SEXP, SEXP);
extern SEXP gRm_item_parameters_ice_fields_from_gamma(SEXP, SEXP);
extern SEXP gRm_gllrm_expected_margins(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"gRm_screen_j_exact_chi_gamma_slices", (DL_FUNC) &gRm_screen_j_exact_chi_gamma_slices, 7},
  {"gRm_screen_j_exact_chi_gamma_trace_slices", (DL_FUNC) &gRm_screen_j_exact_chi_gamma_trace_slices, 7},
  {"gRm_screen_j_exact_chi_slices", (DL_FUNC) &gRm_screen_j_exact_chi_slices, 6},
  {"gRm_screen_j_exact_gamma_slices", (DL_FUNC) &gRm_screen_j_exact_gamma_slices, 4},
  {"gRm_screen_j_conditional_bias_test", (DL_FUNC) &gRm_screen_j_conditional_bias_test, 12},
  {"gRm_screen_j_item_pair_conditional_exact", (DL_FUNC) &gRm_screen_j_item_pair_conditional_exact, 11},
  {"gRm_screen_j_source_random_draws", (DL_FUNC) &gRm_screen_j_source_random_draws, 2},
  /* Compatibility alias retained for older validation callers. */
  {"_gRm_screen_j_source_random_draws", (DL_FUNC) &gRm_screen_j_source_random_draws, 2},
  {"gRm_item_parameters_extended_capabilities", (DL_FUNC) &gRm_item_parameters_extended_capabilities, 0},
  {"gRm_item_parameters_top_ice_field", (DL_FUNC) &gRm_item_parameters_top_ice_field, 2},
  {"gRm_item_parameters_ice_fields_from_gamma", (DL_FUNC) &gRm_item_parameters_ice_fields_from_gamma, 2},
  {"gRm_gllrm_expected_margins", (DL_FUNC) &gRm_gllrm_expected_margins, 4},
  {NULL, NULL, 0}
};

void attribute_visible R_init_gRm(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
