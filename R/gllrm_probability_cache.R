# Implementation-only probability cache.
#
# This file does not implement a separate Pascal algorithm. It memoizes repeated
# GLLRM probability calculations that are otherwise defined by
# source/digram_source_20260817/skunits/skbias22.pas::CalculateBiasedGammaValues2 and the fit loop in
# source/digram_source_20260817/skunits/skbias22.pas::GLLRM_estim. The cache does not change the source algorithm,
# fitted values, update order, or reporting gauge.

#' Internal new gllrm probability cache helper
#'
#' Supports the gllrm probability cache implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param components Internal `components` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gllrm_probability_cache <- function(context, state, components = NULL) {
  list(
    context = context,
    state = state,
    components = components %||% context$ld_components_items %||% gllrm_ld_components(context)$items,
    item = new.env(parent = emptyenv()),
    ld = new.env(parent = emptyenv())
  )
}

#' Internal gllrm probability cache key helper
#'
#' Supports the gllrm probability cache implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_probability_cache_key()` computation result.
#' @keywords internal
#' @noRd
gllrm_probability_cache_key <- function(context, total_score, background_values) {
  dif_backgrounds <- context$dif_background_indices %||% seq_along(background_values)
  if (length(dif_backgrounds) == 0L) {
    return(as.character(total_score))
  }
  paste(c(total_score, background_values[dif_backgrounds]), collapse = "\r")
}

#' Internal gllrm cached item probabilities helper
#'
#' Supports the gllrm probability cache implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param cache Computation cache.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_cached_item_probabilities()` computation result.
#' @keywords internal
#' @noRd
gllrm_cached_item_probabilities <- function(cache, total_score, background_values) {
  key <- gllrm_probability_cache_key(cache$context, total_score, background_values)
  if (!exists(key, envir = cache$item, inherits = FALSE)) {
    assign(
      key,
      gllrm_group_item_probabilities(
        cache$context,
        cache$state,
        total_score = total_score,
        background_values = background_values,
        components = cache$components
      ),
      envir = cache$item
    )
  }
  get(key, envir = cache$item, inherits = FALSE)
}

#' Internal gllrm cached ld probabilities helper
#'
#' Supports the gllrm probability cache implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param cache Computation cache.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_cached_ld_probabilities()` computation result.
#' @keywords internal
#' @noRd
gllrm_cached_ld_probabilities <- function(cache, total_score, background_values) {
  key <- gllrm_probability_cache_key(cache$context, total_score, background_values)
  if (!exists(key, envir = cache$ld, inherits = FALSE)) {
    assign(
      key,
      gllrm_group_ld_probabilities(
        cache$context,
        cache$state,
        total_score = total_score,
        background_values = background_values,
        components = cache$components
      ),
      envir = cache$ld
    )
  }
  get(key, envir = cache$ld, inherits = FALSE)
}
