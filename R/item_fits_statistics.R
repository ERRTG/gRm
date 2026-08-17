#' Internal goodman kruskal gamma helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `goodman_kruskal_gamma()` computation result.
#' @keywords internal
#' @noRd
goodman_kruskal_gamma <- function(tab) {
  source_rc_gamma_stats(tab, include_cells = FALSE)$gamma
}

#' Internal fitted gamma stats helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param expected Internal `expected` value used by this helper.
#' @return The internal `fitted_gamma_stats()` computation result.
#' @keywords internal
#' @noRd
fitted_gamma_stats <- function(expected) {
  # This fitted-gamma variance is the item fit source branch and is
  # intentionally not merged with the exo/screen gamma test variance
  # conventions.
  # Source trace: skbias15.pas::CalculateFittedGAMMA. The source computes
  # gamma = PMQ / PPQ and S1 = 16 / PPQ^4 * sum(tab[i,j] *
  # (Q * AIJ[i,j] - P * DIJ[i,j])^2).
  cd <- gamma_cell_tables(expected)
  ppq <- cd$p + cd$q
  pmq <- cd$p - cd$q
  if (ppq <= 0) {
    return(list(gamma = 0, variance = 0))
  }
  factor <- 16 / (ppq^4)
  variance <- 0
  for (row in seq_len(nrow(expected))) {
    for (col in seq_len(ncol(expected))) {
      m <- cd$q * cd$aij[row, col] - cd$p * cd$dij[row, col]
      variance <- variance + expected[row, col] * m * m
    }
  }
  list(gamma = pmq / ppq, variance = factor * variance)
}

#' Internal item fdr risk helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param p_values Probability values.
#' @return The internal `item_fdr_risk()` computation result.
#' @keywords internal
#' @noRd
item_fdr_risk <- function(p_values) {
  risks <- integer(length(p_values))
  for (alpha_index in seq_along(c(0.05, 0.01, 0.001))) {
    alpha <- c(0.05, 0.01, 0.001)[[alpha_index]]
    critical <- source_bh_critical(p_values, alpha)
    risks[p_values <= critical] <- alpha_index
  }
  risks
}

#' Internal item fit direction helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param rows Rows used by the computation.
#' @return The internal `item_fit_direction()` computation result.
#' @keywords internal
#' @noRd
item_fit_direction <- function(rows) {
  any_flag <- rows$outfit_fdr > 0L | rows$infit_fdr > 0L | rows$gamma_fdr > 0L
  direction <- rep("", nrow(rows))
  direction[any_flag & rows$observed_gamma < rows$expected_gamma] <- "low"
  direction[any_flag & rows$observed_gamma > rows$expected_gamma] <- "high"
  direction
}

#' Internal safe ratio helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param numerator Internal `numerator` value used by this helper.
#' @param denominator Internal `denominator` value used by this helper.
#' @return The internal `safe_ratio()` computation result.
#' @keywords internal
#' @noRd
safe_ratio <- function(numerator, denominator) {
  result <- rep(0, length(numerator))
  ok <- denominator > 0
  result[ok] <- numerator[ok] / denominator[ok]
  result
}

#' Internal safe z helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param observed Internal `observed` value used by this helper.
#' @param expected Internal `expected` value used by this helper.
#' @param variance Internal `variance` value used by this helper.
#' @return The internal `safe_z()` computation result.
#' @keywords internal
#' @noRd
safe_z <- function(observed, expected, variance) {
  z <- rep(0, length(observed))
  ok <- variance > 0
  z[ok] <- (observed[ok] - expected[ok]) / sqrt(variance[ok])
  z
}

#' Internal inf replace helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param x Object or value to process.
#' @param value Value to validate or transform.
#' @return The internal `inf_replace()` computation result.
#' @keywords internal
#' @noRd
inf_replace <- function(x, value) {
  x[!is.finite(x)] <- value
  x
}

#' Internal two sided source normal p helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param z Internal `z` value used by this helper.
#' @return The internal `two_sided_source_normal_p()` computation result.
#' @keywords internal
#' @noRd
two_sided_source_normal_p <- function(z) {
  vapply(z, function(value) pmin(1, 2 * source_tail_norm(abs(value), TRUE)), numeric(1))
}

#' Source fixed-field integer rounding
#'
#' Pascal's fixed-width `:0` numeric formatting rounds half values away from
#' zero for the non-negative item fit frequencies printed by `skbias15.pas`.
#'
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param value Numeric vector.
#' @return Integer-like numeric vector rounded as DIGRAM prints it.
#' @keywords internal
#' @noRd
source_print_round <- function(value) {
  sign(value) * floor(abs(value) + 0.5)
}
