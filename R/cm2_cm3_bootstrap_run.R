#' Run the source-shaped CM2/CM3 parametric bootstrap
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' Generation/finalization trace:
#' `source/digram_source_20260817/skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`
#' and `source/digram_source_20260817/skunits/SKbias8.pas::Finish_random_Gllrm`.
#' @param fit Public fitted model object.
#' @param values Deterministic observed values.
#' @param context Observed-data GLLRM context.
#' @param state Observed fitted state.
#' @param include_three_way Whether to include CM3 margins and aggregate.
#' @param control Normalized bootstrap control.
#' @return The updated values list with a structured `bootstrap` component.
#' @keywords internal
cm2_cm3_apply_parametric_bootstrap <- function(fit,
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

  capability <- cm2_cm3_bootstrap_capability(context, state)
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
      requested_bootstrap_jobs = control$jobs,
      bootstrap_jobs = 0L,
      execution_mode = "not_run",
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

  groups <- cm2_cm3_bootstrap_groups(context)
  distribution_cache <- new_cm2_cm3_bootstrap_distribution_cache(context, state)
  # Source trace: skbias14.CM3_analysis invokes Randomize after
  # Start_random_GLLRM_with_exogene has initialized every probability record
  # and immediately before the bootstrap loop. Initialize the private Delphi
  # stream at the corresponding boundary; an explicit seed is its uint32
  # RandSeed bit pattern and bypasses wall-clock initialization.
  rng <- new_cm2_cm3_bootstrap_rng(control$seed)
  initial_seed <- rng$seed()
  spec <- fit$model %||% fit$spec
  # SKbias8.Estimate_the_GLLRM hard-codes deltalimit = 0.0001 for every
  # bootstrap refit, independently of the observed fit's convergence control.
  tolerance <- 0.0001
  observed_aggregates <- cm2_cm3_aggregate_values(values$tests, context)
  observed_aggregate_rows <- cm2_cm3_bootstrap_aggregate_rows(
    observed_aggregates,
    include_three_way
  )
  observed_margin_rows <- values$tests
  observed_margin_rows$margin_index <- seq_len(nrow(observed_margin_rows))

  # Estimate_GLLRM leaves the loop's stopping-step discrepancy in Pascal's
  # global delta, then CM3_analysis tests that exact value against 0.1. The R
  # fitter separately recomputes state$delta after stopping, so both the source
  # refit path and the no-refit R extension must prefer report_delta here.
  observed_delta <- state$report_delta %||% fit$convergence$report_delta %||%
    state$delta %||% fit$convergence$final_delta %||% 0
  # Start_random_Gllrm_with_exogene.MakeIRTcopy saves the current report-facing
  # IRT structures once. For LD models, reproduce the final source gauge before
  # using that immutable copy as every replicate's starting parameters.
  bootstrap_initial_parameters <- gllrm_output_parameter_state(context, state)
  run_context <- list(
    fit = fit,
    values = values,
    context = context,
    state = state,
    bootstrap_initial_parameters = bootstrap_initial_parameters,
    groups = groups,
    distribution_cache = distribution_cache,
    rng = rng,
    spec = spec,
    tolerance = tolerance,
    observed_delta = observed_delta,
    include_three_way = include_three_way,
    control = control
  )

  # skbias14.CM3_analysis calls Prepare_Extra_M3tests once before DoSomething;
  # each replicate therefore reuses `values$margin_specs` exactly. It does not
  # rerun SelectItems, widen a subset, or change the full model used below.

  # Canonical scd/DIGRAM1f command 197 always calls CM3_analysis with
  # Estimate=true. A FALSE control is the R extension corresponding to the
  # dormant per-sample no-refit branch in DoSomething. It must not be inferred
  # from the Pascal footer: Finish_random_Gllrm restores the saved IRT
  # parameter structures, leaves simulated-data mode, refits the observed
  # data, and writes "Item parameters were reestimated" regardless of the
  # Estimate argument. R never mutates the observed fit and needs no such
  # cleanup refit; `control$reestimate` records only per-sample behavior.

  # Source trace: skbias14.CM3_analysis and
  # SKbias8.Generate_GLLRM_boot_sample_with_exogene generate replicates in one
  # sequential loop using a single mutable RandSeed. Finish that entire phase
  # before any fork so worker scheduling can never change a draw, sample, or
  # later replicate's initial state.
  generated_samples <- cm2_cm3_bootstrap_generate_samples(run_context)
  final_rng_state <- rng$state()
  rng_draws <- rng$draws()

  # The source performs fitting and CM3 analysis immediately after each
  # generated sample. Those routines contain no Random calls, so evaluating
  # the already-generated samples independently is an R execution optimization
  # that preserves every source computation and the replicate result order.
  analysis_context <- run_context
  analysis_context[c("groups", "distribution_cache", "rng")] <- NULL
  bootstrap_jobs <- cm2_cm3_bootstrap_worker_count(
    control$jobs,
    control$nsim
  )
  results <- cm2_cm3_bootstrap_map_samples(
    generated_samples,
    analysis_context,
    bootstrap_jobs
  )
  fit_rows <- lapply(results, `[[`, "fit_row")
  aggregate_rows <- lapply(results, `[[`, "aggregate_rows")
  margin_rows <- lapply(results, `[[`, "margin_rows")
  samples <- if (isTRUE(control$keep_samples)) generated_samples else list()

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
  aggregate_summary <- cm2_cm3_bootstrap_summarize_rows(
    observed_aggregate_rows,
    aggregate_table,
    key = c("diagnostic", "kind", "background_label"),
    nused = nused
  )
  margin_summary <- cm2_cm3_bootstrap_summarize_rows(
    observed_margin_rows,
    margin_table,
    key = "margin_index",
    nused = nused
  )
  values <- cm2_cm3_attach_bootstrap_summaries(
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
    final_rng_state = final_rng_state,
    rng_draws = rng_draws,
    rng = "delphi4_lcg_134775813_uint32",
    generation_mode = unique(fit_table$generation_mode),
    reestimate = control$reestimate,
    requested_bootstrap_jobs = control$jobs,
    bootstrap_jobs = bootstrap_jobs,
    execution_mode = if (bootstrap_jobs > 1L) {
      "ordered_fork_parallel"
    } else {
      "serial"
    },
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

#' Generate all source-shaped CM2/CM3 bootstrap samples in RNG order
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`
#' and `source/digram_source_20260817/skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param run_context Bootstrap inputs including the single mutable Delphi RNG.
#' @return A replicate-ordered list of generated samples.
#' @keywords internal
#' @noRd
cm2_cm3_bootstrap_generate_samples <- function(run_context) {
  samples <- vector("list", run_context$control$nsim)
  for (replicate in seq_len(run_context$control$nsim)) {
    samples[[replicate]] <- cm2_cm3_bootstrap_generate_sample(
      run_context$context,
      run_context$groups,
      run_context$distribution_cache,
      run_context$rng
    )
  }
  samples
}

#' Resolve the effective bootstrap analysis worker count
#'
#' Forked workers are available only on POSIX platforms. Generation is always
#' serial, and a request cannot use more workers than generated replicates or
#' detected logical cores.
#' @param jobs Normalized requested worker count.
#' @param nsim Number of generated bootstrap samples.
#' @return A positive integer worker count.
#' @keywords internal
#' @noRd
cm2_cm3_bootstrap_worker_count <- function(jobs, nsim) {
  if (.Platform$OS.type != "unix") {
    return(1L)
  }
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) {
    detected <- 1L
  }
  as.integer(max(1L, min(jobs, nsim, detected)))
}

#' Analyze generated CM2/CM3 bootstrap samples in stable order
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`
#' and `source/digram_source_20260817/skunits/SKbias8.pas::Estimate_the_GLLRM`.
#' The source computations in this phase are RNG-free. On POSIX, ordered fork
#' mapping changes only their scheduling; `mc.set.seed = FALSE` also prevents
#' the R wrapper from creating worker RNG substreams that the source never had.
#' @param samples Replicate-ordered generated samples.
#' @param run_context Immutable fitting and analysis inputs without the RNG.
#' @param jobs Effective worker count.
#' @return A replicate-ordered list of analyzed results.
#' @keywords internal
#' @noRd
cm2_cm3_bootstrap_map_samples <- function(samples, run_context, jobs) {
  indices <- seq_along(samples)
  analyze_one <- function(replicate) {
    tryCatch(
      cm2_cm3_bootstrap_replicate(
        replicate,
        samples[[replicate]],
        run_context
      ),
      error = function(error) {
        structure(
          list(replicate = replicate, message = conditionMessage(error)),
          class = "gRm_bootstrap_worker_error"
        )
      }
    )
  }
  results <- if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(
      indices,
      analyze_one,
      mc.cores = jobs,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(indices, analyze_one)
  }
  failed <- which(vapply(
    results,
    inherits,
    logical(1L),
    what = "gRm_bootstrap_worker_error"
  ))
  if (length(failed)) {
    failure <- results[[failed[[1L]]]]
    stop(
      "CM2/CM3 bootstrap replicate ", failure$replicate,
      " failed: ", failure$message,
      call. = FALSE
    )
  }
  results
}

#' Analyze one source-shaped CM2/CM3 bootstrap replicate
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias8.pas::Estimate_the_GLLRM`.
#' Analysis trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' Mathematical step: optionally perform one fresh GLLRM fit on an already
#' generated conditional sample, apply the strict source stopping-delta
#' acceptance rule, and evaluate all requested diagnostic margins without
#' changing their order.
#' @param replicate One-based bootstrap replicate index.
#' @param sample One sequentially generated sample.
#' @param run_context Immutable bootstrap fitting and analysis inputs.
#' @return A list containing fit, aggregate, and margin rows.
#' @keywords internal
#' @noRd
cm2_cm3_bootstrap_replicate <- function(replicate, sample, run_context) {
  control <- run_context$control
  sample_bundle <- cm2_cm3_bootstrap_bundle(run_context$fit, sample)
  sample_context <- build_gllrm_context(run_context$spec, sample_bundle)
  sample_state <- run_context$state
  fit_error <- NA_character_

  # This is CM3_analysis.DoSomething's `If Estimate then` branch. Canonical
  # command 197 passes true; FALSE is an explicit package extension in which
  # CM3_tests still uses the supplied observed parameters. The unconditional
  # observed-data refit/footer in Finish_random_Gllrm is later cleanup, not a
  # simulated-sample refit.
  if (isTRUE(control$reestimate)) {
    fitted <- tryCatch(
      fit_gllrm(
        run_context$spec,
        max_step = control$max_step,
        max_delta = run_context$tolerance,
        bundle = sample_bundle,
        initial_parameters = run_context$bootstrap_initial_parameters
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
    # Source trace: SKbias8.Estimate_GLLRM/Estimate_the_GLLRM returns with the
    # stopping-step global delta; skbias14.CM3_analysis immediately accepts on
    # that value. state$delta is a later R-only recomputation and is not used.
    as.numeric(sample_state$report_delta %||% NA_real_)
  } else {
    as.numeric(run_context$observed_delta)
  }
  accepted <- !is.null(sample_state) && is.finite(final_delta) &&
    final_delta < control$acceptance_delta
  analysis <- if (is.null(sample_state)) {
    NULL
  } else {
    tryCatch(
      cm2_cm3_bootstrap_analyze(
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
    aggregate_rows <- cm2_cm3_bootstrap_aggregate_rows(
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
    margin_rows = margin_rows
  )
}
