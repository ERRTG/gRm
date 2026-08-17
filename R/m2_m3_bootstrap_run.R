#' Run the source-shaped CM2/CM3 parametric bootstrap
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param fit Public fitted model object.
#' @param values Deterministic observed values.
#' @param context Observed-data GLLRM context.
#' @param state Observed fitted state.
#' @param include_three_way Whether to include M3 margins and aggregate.
#' @param control Normalized bootstrap control.
#' @return The updated values list with a structured `bootstrap` component.
#' @keywords internal
m2_m3_apply_parametric_bootstrap <- function(fit,
                                              values,
                                              context,
                                              state,
                                              include_three_way,
                                              control) {
  if (!isTRUE(control$enabled)) {
    values$bootstrap <- list(
      enabled = FALSE,
      possible = NA,
      source_status = "not_requested",
      nsim = 0L,
      nused = 0L
    )
    return(values)
  }

  capability <- m2_m3_bootstrap_capability(context)
  if (!isTRUE(capability$possible)) {
    values$bootstrap <- list(
      enabled = TRUE,
      possible = FALSE,
      source_status = "bootstrap_not_possible",
      reasons = capability$reasons,
      capability = capability,
      nsim = control$nsim,
      nused = 0L,
      seed = control$seed,
      reestimate = control$reestimate,
      resample_score_distribution = control$resample_score_distribution,
      score_distribution_mode = "fixed_observed_source_noop_branch",
      replicates = data.frame(),
      aggregate_replicates = data.frame(),
      margin_replicates = data.frame(),
      aggregate_summary = data.frame(),
      margin_summary = data.frame(),
      samples = list()
    )
    return(values)
  }

  groups <- m2_m3_bootstrap_groups(context)
  rng <- new_m2_m3_bootstrap_rng(control$seed)
  initial_seed <- rng$seed()
  distribution_cache <- new_m2_m3_bootstrap_distribution_cache(context, state)
  spec <- fit$model %||% fit$spec
  tolerance <- fit$convergence$tolerance %||% fit$convergence$max_delta %||% 0.0001
  observed_aggregates <- m2_m3_aggregate_values(values$tests)
  observed_aggregate_rows <- m2_m3_bootstrap_aggregate_rows(
    observed_aggregates,
    include_three_way
  )
  observed_margin_rows <- values$tests
  observed_margin_rows$margin_index <- seq_len(nrow(observed_margin_rows))

  fit_rows <- vector("list", control$nsim)
  aggregate_rows <- vector("list", control$nsim)
  margin_rows <- vector("list", control$nsim)
  samples <- if (isTRUE(control$keep_samples)) vector("list", control$nsim) else list()

  observed_delta <- state$delta %||% fit$convergence$final_delta %||%
    fit$convergence$report_delta %||% 0
  run_context <- list(
    fit = fit,
    values = values,
    context = context,
    state = state,
    groups = groups,
    distribution_cache = distribution_cache,
    rng = rng,
    spec = spec,
    tolerance = tolerance,
    observed_delta = observed_delta,
    include_three_way = include_three_way,
    control = control
  )
  for (replicate in seq_len(control$nsim)) {
    result <- m2_m3_bootstrap_replicate(replicate, run_context)
    fit_rows[[replicate]] <- result$fit_row
    aggregate_rows[[replicate]] <- result$aggregate_rows
    margin_rows[[replicate]] <- result$margin_rows
    if (isTRUE(control$keep_samples)) {
      samples[[replicate]] <- result$sample
    }
  }

  fit_table <- do.call(rbind, fit_rows)
  aggregate_rows <- Filter(is.data.frame, aggregate_rows)
  margin_rows <- Filter(is.data.frame, margin_rows)
  aggregate_table <- if (length(aggregate_rows)) {
    do.call(rbind_fill, aggregate_rows)
  } else {
    data.frame()
  }
  margin_table <- if (length(margin_rows)) {
    do.call(rbind_fill, margin_rows)
  } else {
    data.frame()
  }
  nused <- sum(fit_table$accepted)
  aggregate_summary <- m2_m3_bootstrap_summarize_rows(
    observed_aggregate_rows,
    aggregate_table,
    key = c("diagnostic", "kind", "background_label"),
    nused = nused
  )
  margin_summary <- m2_m3_bootstrap_summarize_rows(
    observed_margin_rows,
    margin_table,
    key = "margin_index",
    nused = nused
  )
  values <- m2_m3_attach_bootstrap_summaries(
    values,
    aggregate_summary,
    margin_summary,
    include_three_way
  )
  values$bootstrap <- list(
    enabled = TRUE,
    possible = TRUE,
    source_status = "source_shaped_parametric_bootstrap",
    nsim = control$nsim,
    nused = as.integer(nused),
    seed = initial_seed,
    final_rng_state = rng$state(),
    rng_draws = rng$draws(),
    rng = "park_miller_48271_private_validation_stream",
    generation_mode = unique(fit_table$generation_mode),
    reestimate = control$reestimate,
    max_step = control$max_step,
    max_delta = tolerance,
    acceptance_delta = control$acceptance_delta,
    exceedance_direction = "simulated_p_value <= observed_p_value",
    resample_score_distribution = control$resample_score_distribution,
    score_distribution_mode = "fixed_observed_source_noop_branch",
    capability = capability,
    score_exogenous_groups = groups,
    replicates = fit_table,
    aggregate_replicates = aggregate_table,
    margin_replicates = margin_table,
    aggregate_summary = aggregate_summary,
    margin_summary = margin_summary,
    samples = samples
  )
  values
}

#' Run one source-shaped CM2/CM3 bootstrap replicate
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' Mathematical step: generate one conditional sample, optionally perform one
#' fresh GLLRM fit, apply the strict final-delta acceptance rule, and evaluate
#' all requested diagnostic margins without changing their traversal order.
#' @param replicate One-based bootstrap replicate index.
#' @param run_context Immutable bootstrap inputs plus the private mutable RNG.
#' @return A list containing fit, aggregate, margin, and generated-sample rows.
#' @keywords internal
#' @noRd
m2_m3_bootstrap_replicate <- function(replicate, run_context) {
  control <- run_context$control
  sample <- m2_m3_bootstrap_generate_sample(
    run_context$context,
    run_context$groups,
    run_context$distribution_cache,
    run_context$rng
  )
  sample_bundle <- m2_m3_bootstrap_bundle(run_context$fit, sample)
  sample_context <- build_gllrm_context(run_context$spec, sample_bundle)
  sample_state <- run_context$state
  fit_error <- NA_character_

  if (isTRUE(control$reestimate)) {
    fitted <- tryCatch(
      fit_gllrm(
        run_context$spec,
        max_step = control$max_step,
        max_delta = run_context$tolerance,
        bundle = sample_bundle
      ),
      error = function(error) error
    )
    if (inherits(fitted, "error")) {
      fit_error <- conditionMessage(fitted)
      sample_state <- NULL
    } else {
      sample_context <- fitted$context
      sample_state <- fitted$state
    }
  }

  final_delta <- if (is.null(sample_state)) {
    NA_real_
  } else if (isTRUE(control$reestimate)) {
    as.numeric(sample_state$delta %||% NA_real_)
  } else {
    as.numeric(run_context$observed_delta)
  }
  accepted <- !is.null(sample_state) && is.finite(final_delta) &&
    final_delta < control$acceptance_delta
  analysis <- if (is.null(sample_state)) {
    NULL
  } else {
    tryCatch(
      m2_m3_bootstrap_analyze(
        sample_context,
        sample_state,
        specs = run_context$values$margin_specs,
        score_cuts = run_context$values$score_cuts
      ),
      error = function(error) error
    )
  }
  if (inherits(analysis, "error")) {
    fit_error <- conditionMessage(analysis)
    analysis <- NULL
    accepted <- FALSE
  }

  fit_row <- data.frame(
    replicate = replicate,
    accepted = accepted,
    reestimated = control$reestimate,
    iterations = if (is.null(sample_state)) NA_integer_ else as.integer(sample_state$n_step %||% NA_integer_),
    report_delta = if (is.null(sample_state)) NA_real_ else as.numeric(sample_state$report_delta %||% NA_real_),
    final_delta = final_delta,
    converged = if (is.null(sample_state)) NA else isTRUE(sample_state$converged),
    stop_reason = if (is.null(sample_state)) NA_character_ else as.character(sample_state$stop_reason %||% NA_character_),
    error = fit_error,
    n_records = nrow(sample$items),
    rng_start_state = sample$start_state,
    rng_final_state = sample$final_state,
    rng_draws = sample$draws,
    generation_mode = sample$generation_mode,
    stringsAsFactors = FALSE
  )
  aggregate_rows <- margin_rows <- NULL
  if (!is.null(analysis)) {
    aggregate_rows <- m2_m3_bootstrap_aggregate_rows(
      analysis$aggregates,
      run_context$include_three_way
    )
    aggregate_rows$replicate <- replicate
    aggregate_rows$accepted <- accepted

    margin_rows <- analysis$tests
    margin_rows$margin_index <- seq_len(nrow(margin_rows))
    margin_rows$replicate <- replicate
    margin_rows$accepted <- accepted
  }

  list(
    fit_row = fit_row,
    aggregate_rows = aggregate_rows,
    margin_rows = margin_rows,
    sample = sample
  )
}
