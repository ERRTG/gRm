#' Evaluate deterministic CM2/CM3 values in one generated sample
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Generated-sample GLLRM context.
#' @param state Original or re-estimated GLLRM state.
#' @param specs Prepared observed-data margin specifications.
#' @param score_cuts Resolved diagnostic score cuts.
#' @return A list containing per-margin and aggregate statistics.
#' @keywords internal
m2_m3_bootstrap_analyze <- function(context, state, specs, score_cuts) {
  score_group_lookup <- m2_m3_score_group_lookup(context, score_cuts)
  probability_cache <- new_m2_m3_focus_probability_cache(context, state)
  tests <- if (length(specs)) {
    do.call(rbind, lapply(specs, function(spec) {
      m2_m3_analyze_margin(context, state, spec, score_group_lookup, probability_cache)
    }))
  } else {
    data.frame()
  }
  rownames(tests) <- NULL
  list(tests = tests, aggregates = m2_m3_aggregate_values(tests))
}

#' Build one long aggregate-statistic table
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param aggregates Aggregate value list from `m2_m3_aggregate_values()`.
#' @param include_three_way Whether the public result includes M3.
#' @return A data frame containing M2, optional M3, item-trait, and invariance.
#' @keywords internal
m2_m3_bootstrap_aggregate_rows <- function(aggregates, include_three_way) {
  rows <- list(aggregates$m2)
  if (isTRUE(include_three_way)) {
    rows <- c(rows, list(aggregates$m3))
  }
  rows <- c(rows, list(aggregates$item_trait))
  out <- do.call(rbind, rows)
  out$kind <- "aggregate"
  out$background_label <- NA_character_
  out$background_name <- NA_character_
  invariance <- aggregates$invariance
  if (is.data.frame(invariance) && nrow(invariance)) {
    extra <- data.frame(
      diagnostic = paste0("Invariance: ", invariance$background_name),
      chi_square = invariance$chi_square,
      degrees_of_freedom = invariance$degrees_of_freedom,
      p_value = invariance$p_value,
      kind = "invariance",
      background_label = invariance$background_label,
      background_name = invariance$background_name,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, extra)
  }
  rownames(out) <- NULL
  out
}

#' Summarize source-direction bootstrap exceedance counts
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param observed Observed statistic rows.
#' @param simulated Accepted and rejected simulated statistic rows.
#' @param key Key columns identifying corresponding rows.
#' @param nused Number of accepted samples.
#' @return Observed rows with extreme count, denominator, and bootstrap p-value.
#' @keywords internal
m2_m3_bootstrap_summarize_rows <- function(observed, simulated, key, nused) {
  out <- observed
  out$bootstrap_extreme_count <- 0L
  out$bootstrap_nused <- as.integer(nused)
  out$bootstrap_p_value <- NA_real_
  if (!nrow(out) || nused <= 0L) {
    return(out)
  }
  for (row in seq_len(nrow(out))) {
    keep <- simulated$accepted
    for (column in key) {
      value <- out[[column]][[row]]
      keep <- keep & if (is.na(value)) {
        is.na(simulated[[column]])
      } else {
        !is.na(simulated[[column]]) & simulated[[column]] == value
      }
    }
    extreme <- sum(keep & simulated$p_value <= out$p_value[[row]], na.rm = TRUE)
    out$bootstrap_extreme_count[[row]] <- as.integer(extreme)
    out$bootstrap_p_value[[row]] <- extreme / nused
  }
  out
}

#' Attach bootstrap calibration columns to deterministic result tables
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param values Deterministic M2/M3 values.
#' @param aggregate_summary Aggregate bootstrap summary.
#' @param margin_summary Per-margin bootstrap summary.
#' @param include_three_way Whether `values` contains an M3 aggregate.
#' @return Updated values list.
#' @keywords internal
m2_m3_attach_bootstrap_summaries <- function(values,
                                              aggregate_summary,
                                              margin_summary,
                                              include_three_way) {
  attach_aggregate <- function(table, diagnostic) {
    match_row <- match(diagnostic, aggregate_summary$diagnostic)
    table$bootstrap_extreme_count <- aggregate_summary$bootstrap_extreme_count[[match_row]]
    table$bootstrap_nused <- aggregate_summary$bootstrap_nused[[match_row]]
    table$bootstrap_p_value <- aggregate_summary$bootstrap_p_value[[match_row]]
    table
  }
  if (isTRUE(include_three_way)) {
    values$m2 <- attach_aggregate(values$m2, "M2")
    values$m3 <- attach_aggregate(values$m3, "M3")
  } else {
    values$aggregate <- attach_aggregate(values$aggregate, "M2")
  }
  values$item_trait <- attach_aggregate(values$item_trait, "Item-trait")
  if (nrow(values$invariance)) {
    matches <- match(values$invariance$background_label, aggregate_summary$background_label)
    values$invariance$bootstrap_extreme_count <- aggregate_summary$bootstrap_extreme_count[matches]
    values$invariance$bootstrap_nused <- aggregate_summary$bootstrap_nused[matches]
    values$invariance$bootstrap_p_value <- aggregate_summary$bootstrap_p_value[matches]
  } else {
    values$invariance$bootstrap_extreme_count <- integer()
    values$invariance$bootstrap_nused <- integer()
    values$invariance$bootstrap_p_value <- numeric()
  }
  if (nrow(values$tests)) {
    values$tests$bootstrap_extreme_count <- margin_summary$bootstrap_extreme_count
    values$tests$bootstrap_nused <- margin_summary$bootstrap_nused
    values$tests$bootstrap_p_value <- margin_summary$bootstrap_p_value
  }
  values
}
