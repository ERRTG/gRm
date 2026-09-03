#' Evaluate deterministic CM2/CM3 values in one generated sample
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_tests`,
#' called from the nested `CM3_analysis.DoSomething` bootstrap step.
#' @param context Generated-sample GLLRM context.
#' @param state Original or re-estimated GLLRM state.
#' @param specs Margin specifications prepared once from the observed model and
#'   reused for every generated sample, as in `CM3_analysis`.
#' @param score_cuts Resolved diagnostic score cuts.
#' @return A list containing per-margin and aggregate statistics.
#' @keywords internal
cm2_cm3_bootstrap_analyze <- function(context, state, specs, score_cuts) {
  # CM3_analysis prepares its M3tests arrays before the bootstrap loop. Each
  # DoSomething call then invokes CM3_tests on SimIRTdata using that unchanged
  # margin order; rebuild only sample-dependent counts and probabilities here.
  score_group_lookup <- cm2_cm3_score_group_lookup(context, score_cuts)
  probability_cache <- new_cm2_cm3_focus_probability_cache(context, state)
  tests <- if (length(specs)) {
    do.call(rbind, lapply(specs, function(spec) {
      cm2_cm3_analyze_margin(context, state, spec, score_group_lookup, probability_cache)
    }))
  } else {
    data.frame()
  }
  rownames(tests) <- NULL
  # CM3_tests derives CM2, CM3, item-trait, and per-exogenous invariance from
  # the same ordered margin rows. Keep that common derivation for every sample.
  list(tests = tests, aggregates = cm2_cm3_aggregate_values(tests, context))
}

#' Build one long aggregate-statistic table
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' This is an R storage form for the aggregate values whose comparisons update
#' scalar `Nlarger2`, `Nlarger`, `NlargerIT`, and vector
#' `NlargerInvariance` counters in `DoSomething` before `Summarize` emits them.
#' @param aggregates Aggregate value list from `cm2_cm3_aggregate_values()`.
#' @param include_three_way Whether the public result includes CM3.
#' @return A data frame containing CM2, optional CM3, item-trait, and invariance.
#' @keywords internal
cm2_cm3_bootstrap_aggregate_rows <- function(aggregates, include_three_way) {
  rows <- list(aggregates$cm2)
  if (isTRUE(include_three_way)) {
    rows <- c(rows, list(aggregates$cm3))
  }
  rows <- c(rows, list(aggregates$item_trait))
  # Pascal owns separate scalar values and counters for CM2, CM3, and
  # item-trait. The R long table changes only storage shape so one source
  # comparison rule can be used.
  out <- do.call(rbind, rows)
  out$kind <- "aggregate"
  out$background_label <- NA_character_
  out$background_name <- NA_character_
  invariance <- aggregates$invariance
  if (is.data.frame(invariance) && nrow(invariance)) {
    # CM3_analysis loops over exogenous variables in source order when updating
    # and printing NlargerInvariance; preserve that order in the appended rows.
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
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' Its nested `DoSomething` increments counters only when `delta < 0.1` and
#' when a simulated asymptotic p-value is no larger than its observed value;
#' nested `Summarize` divides every counter by `nused` without a plus-one
#' correction.
#' @param observed Observed statistic rows.
#' @param simulated Accepted and rejected simulated statistic rows.
#' @param key Key columns identifying corresponding rows.
#' @param nused Number of accepted samples.
#' @return Observed rows with extreme count, denominator, and bootstrap p-value.
#' @keywords internal
cm2_cm3_bootstrap_summarize_rows <- function(observed, simulated, key, nused) {
  out <- observed
  # CM3_analysis initializes Nlarger/Nlarger2/NlargerIT, the invariance vector,
  # and every M3tests0^.antal counter to zero before entering its sample loop.
  out$bootstrap_extreme_count <- 0L
  out$bootstrap_nused <- as.integer(nused)
  out$bootstrap_p_value <- NA_real_
  if (!nrow(out) || nused <= 0L) {
    return(out)
  }
  for (row in seq_len(nrow(out))) {
    # Pascal's acceptance guard encloses all counter increments. `accepted`
    # records the exact strict global-delta test performed in DoSomething.
    keep <- simulated$accepted
    for (column in key) {
      value <- out[[column]][[row]]
      keep <- keep & if (is.na(value)) {
        is.na(simulated[[column]])
      } else {
        !is.na(simulated[[column]]) & simulated[[column]] == value
      }
    }
    # Source direction is p_sim <= p_observed. Do not reverse the comparison,
    # compare chi-square values, or add a Monte Carlo continuity correction.
    extreme <- sum(keep & simulated$p_value <= out$p_value[[row]], na.rm = TRUE)
    out$bootstrap_extreme_count[[row]] <- as.integer(extreme)
    out$bootstrap_p_value[[row]] <- extreme / nused
  }
  out
}

#' Attach bootstrap calibration columns to deterministic result tables
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' Pascal's nested `Summarize` prints the same counters in margin, CM2, CM3,
#' item-trait, and invariance sections; this helper attaches them to the
#' corresponding structured R tables without applying its display-only filter.
#' @param values Deterministic CM2/CM3 values.
#' @param aggregate_summary Aggregate bootstrap summary.
#' @param margin_summary Per-margin bootstrap summary.
#' @param include_three_way Whether `values` contains an CM3 aggregate.
#' @return Updated values list.
#' @keywords internal
cm2_cm3_attach_bootstrap_summaries <- function(values,
                                              aggregate_summary,
                                              margin_summary,
                                              include_three_way) {
  attach_aggregate <- function(table, diagnostic) {
    # These labels replace Pascal's separate scalar variables; matching them
    # changes representation only and retains the source aggregate identities.
    match_row <- match(diagnostic, aggregate_summary$diagnostic)
    table$bootstrap_extreme_count <- aggregate_summary$bootstrap_extreme_count[[match_row]]
    table$bootstrap_nused <- aggregate_summary$bootstrap_nused[[match_row]]
    table$bootstrap_p_value <- aggregate_summary$bootstrap_p_value[[match_row]]
    table
  }
  if (isTRUE(include_three_way)) {
    values$cm2 <- attach_aggregate(values$cm2, "CM2")
    values$cm3 <- attach_aggregate(values$cm3, "CM3")
  } else {
    values$aggregate <- attach_aggregate(values$aggregate, "CM2")
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
    # CM3_analysis.Summarize may print only rows with bootstrap p <= 0.05 when
    # PrintAll is false. Retain every prepared row programmatically; rendering
    # or oracle extraction owns that source display policy.
    values$tests$bootstrap_extreme_count <- margin_summary$bootstrap_extreme_count
    values$tests$bootstrap_nused <- margin_summary$bootstrap_nused
    values$tests$bootstrap_p_value <- margin_summary$bootstrap_p_value
  }
  values
}
