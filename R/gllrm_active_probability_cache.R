new_active_gllrm_probability_cache <- function(context, state, components = NULL) {
  list(
    context = context,
    state = state,
    components = components %||% context$ld_components_items %||% gllrm_ld_components(context)$items,
    item = new.env(parent = emptyenv()),
    ld = new.env(parent = emptyenv())
  )
}

active_gllrm_probability_cache_key <- function(context, total_score, background_values) {
  active <- context$active_background_indices %||% seq_along(background_values)
  if (length(active) == 0L) {
    return(as.character(total_score))
  }
  paste(c(total_score, background_values[active]), collapse = "\r")
}

active_gllrm_cached_item_probabilities <- function(cache, total_score, background_values) {
  key <- active_gllrm_probability_cache_key(cache$context, total_score, background_values)
  if (!exists(key, envir = cache$item, inherits = FALSE)) {
    assign(
      key,
      active_gllrm_group_item_probabilities(
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

active_gllrm_cached_ld_probabilities <- function(cache, total_score, background_values) {
  key <- active_gllrm_probability_cache_key(cache$context, total_score, background_values)
  if (!exists(key, envir = cache$ld, inherits = FALSE)) {
    assign(
      key,
      active_gllrm_group_ld_probabilities(
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
