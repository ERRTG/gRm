# Implementation-only probability cache.
#
# This file does not implement a separate Pascal algorithm. It memoizes repeated
# GLLRM probability calculations that are otherwise defined by
# source/GLLRM_ESTIM.txt::CalculateBiasedGammaValues2 and the fit loop in
# source/GLLRM_ESTIM.txt::GLLRM_estim. The cache does not change the source algorithm,
# fitted values, update order, or reporting gauge.

new_gllrm_probability_cache <- function(context, state, components = NULL) {
  list(
    context = context,
    state = state,
    components = components %||% context$ld_components_items %||% gllrm_ld_components(context)$items,
    item = new.env(parent = emptyenv()),
    ld = new.env(parent = emptyenv())
  )
}

gllrm_probability_cache_key <- function(context, total_score, background_values) {
  dif_backgrounds <- context$dif_background_indices %||% seq_along(background_values)
  if (length(dif_backgrounds) == 0L) {
    return(as.character(total_score))
  }
  paste(c(total_score, background_values[dif_backgrounds]), collapse = "\r")
}

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
