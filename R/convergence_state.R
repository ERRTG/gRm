#' Construct the stable public convergence state
#'
#' `report_delta` is the source discrepancy retained at the stopping step;
#' `final_delta` is the diagnostic discrepancy after the final expected-margin
#' recomputation. Both are sufficient-count discrepancies, not parameter
#' movements. The attempted block identifies the numerical attempt whose
#' values were used in the report.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param converged Final source convergence flag.
#' @param iterations Number of completed IPF steps.
#' @param report_delta Source report/stopping discrepancy.
#' @param final_delta Final recomputed diagnostic discrepancy.
#' @param tolerance Requested strict source discrepancy tolerance.
#' @param max_step Requested maximum number of steps.
#' @param stop_reason Stable source stop-reason label.
#' @param source_converged_before_post_acceptance Source flag before any
#'   documented post-stop acceptance rule.
#' @param post_stop_accepted Whether the source post-stop rule changed a false
#'   flag to true.
#' @param report_value_source Label identifying the attempt used for output.
#' @return A `gRm_convergence_state` list.
#' @keywords internal
new_gRm_convergence_state <- function(converged,
                                      iterations,
                                      report_delta,
                                      final_delta,
                                      tolerance,
                                      max_step,
                                      stop_reason,
                                      source_converged_before_post_acceptance,
                                      post_stop_accepted,
                                      report_value_source = "attempted_fit") {
  converged <- isTRUE(converged)
  attempted <- list(
    iterations = as.integer(iterations),
    report_delta = as.numeric(report_delta),
    final_delta = as.numeric(final_delta),
    converged = converged,
    stop_reason = as.character(stop_reason)
  )
  out <- list(
    schema = "gRm-convergence-state-v1",
    converged = converged,
    source_converged = converged,
    source_converged_before_post_acceptance =
      isTRUE(source_converged_before_post_acceptance),
    iterations = as.integer(iterations),
    report_delta = as.numeric(report_delta),
    final_delta = as.numeric(final_delta),
    tolerance = as.numeric(tolerance),
    max_step = as.integer(max_step),
    stop_reason = as.character(stop_reason),
    post_stop_accepted = isTRUE(post_stop_accepted),
    attempted = attempted,
    report_value_source = as.character(report_value_source),
    # Transition aliases retained for code written against gRm 0.x.
    delta = as.numeric(report_delta),
    max_delta = as.numeric(tolerance)
  )
  class(out) <- c("gRm_convergence_state", "list")
  out
}

#' Convert a base Rasch fit to the public convergence schema
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param fit Low-level result from `fit_rasch_base()`.
#' @param max_step Requested maximum number of steps.
#' @param max_delta Requested strict discrepancy tolerance.
#' @return A `gRm_convergence_state` list.
#' @keywords internal
gRm_rasch_convergence_state <- function(fit, max_step, max_delta) {
  stop_reason <- if ((fit$report_delta %||% Inf) < max_delta) {
    "delta_below_tolerance"
  } else {
    "max_step"
  }
  new_gRm_convergence_state(
    converged = fit$converged,
    iterations = fit$n_step,
    report_delta = fit$report_delta,
    final_delta = fit$delta,
    tolerance = max_delta,
    max_step = max_step,
    stop_reason = stop_reason,
    source_converged_before_post_acceptance = fit$converged,
    post_stop_accepted = FALSE
  )
}

#' Convert a GLLRM state to the public convergence schema
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param state Low-level state returned by `fit_gllrm()`.
#' @param max_step Requested maximum number of steps.
#' @param max_delta Requested strict discrepancy tolerance.
#' @return A `gRm_convergence_state` list.
#' @keywords internal
gRm_gllrm_convergence_state <- function(state, max_step, max_delta) {
  before <- isTRUE(state$convergence_before_final_acceptance)
  final <- isTRUE(state$converged)
  new_gRm_convergence_state(
    converged = final,
    iterations = state$n_step,
    report_delta = state$report_delta,
    final_delta = state$delta,
    tolerance = max_delta,
    max_step = max_step,
    stop_reason = state$stop_reason,
    source_converged_before_post_acceptance = before,
    post_stop_accepted = !before && final
  )
}
