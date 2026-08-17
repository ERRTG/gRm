#' Initialize DIGRAM GLLRM stopping control
#'
#' Creates the scalar/history state used by the source IPF stopping machine.
#' The initial discrepancy is `Nvalid` and history indices -4 through 0 are
#' seeded with 9999 exactly as in `Estimate_LL_parameters`.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param n_valid Number of estimation records.
#' @param max_step Source maximum IPF step.
#' @param max_delta Source strict convergence tolerance.
#' @return An internal GLLRM stopping-control list.
#' @keywords internal
source_gllrm_control_state <- function(n_valid, max_step, max_delta) {
  list(
    n_step = 0L,
    delta = as.numeric(n_valid),
    previous_delta = NA_real_,
    delta_history = stats::setNames(rep(9999, 5L), as.character(-4:0)),
    finish = FALSE,
    n_finish = 0L,
    convergence = TRUE,
    recurring = FALSE,
    stop_reason = NA_character_,
    initial_delta = NA_real_,
    min_delta = NA_real_,
    min_delta_step = NA_integer_,
    max_step = as.integer(max_step),
    max_delta = as.numeric(max_delta)
  )
}

#' Read one indexed source delta-history value
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param control A source GLLRM stopping-control list.
#' @param step Integer source history index.
#' @return The stored numeric discrepancy.
#' @keywords internal
source_gllrm_history_value <- function(control, step) {
  value <- unname(control$delta_history[as.character(as.integer(step))])
  if (length(value) != 1L || is.na(value)) {
    stop("Missing source GLLRM delta-history value for step ", step, ".", call. = FALSE)
  }
  value
}

#' Test the source recurring-delta condition
#'
#' Implements `source/GLLRM_ESTIM.txt::RecurringDeltaValues` with its exact
#' `Nstep > 5`, `>= 0.00001`, and current-delta `>= 2` boundaries.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param control A source GLLRM stopping-control list after recording a step.
#' @return A single logical value.
#' @keywords internal
source_gllrm_recurring_delta_values <- function(control) {
  if (control$n_step <= 5L) {
    return(FALSE)
  }
  current <- source_gllrm_history_value(control, control$n_step)
  two_back <- source_gllrm_history_value(control, control$n_step - 2L)
  four_back <- source_gllrm_history_value(control, control$n_step - 4L)
  (current - two_back >= 0.00001) &&
    (current - four_back >= 0.00001) &&
    (current >= 2)
}

#' Apply the exact source GLLRM stop decision
#'
#' The recurrence guard precedes `Iteration_stop`; the remaining branches keep
#' the source expression order and strict comparisons. The returned convergence
#' flag is the pre-finalization source flag.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param control A source GLLRM stopping-control list after recording a step.
#' @return A list with `stop`, `reason`, `convergence`, and `recurring`.
#' @keywords internal
source_gllrm_stop_decision <- function(control) {
  recurring <- source_gllrm_recurring_delta_values(control)
  convergence <- isTRUE(control$convergence)
  stop <- FALSE
  reason <- NA_character_

  # Source trace: source/GLLRM_ESTIM.txt::Estimate_LL_parameters checks this
  # recurrence branch before calling Iteration_stop and only after step 50.
  if (recurring && control$n_step > 50L) {
    stop <- TRUE
    reason <- "recurring_delta_values"
    if (control$delta > control$max_delta) {
      convergence <- FALSE
    }
  } else if (control$delta < control$max_delta) {
    stop <- TRUE
    reason <- "delta_below_tolerance"
  } else if (control$n_step == control$max_step) {
    stop <- TRUE
    reason <- "max_step"
    if (control$delta > control$max_delta) {
      convergence <- FALSE
    }
  } else if (control$n_step %% 1000L == 0L) {
    stop <- TRUE
    reason <- "periodic_1000"
    if (control$delta > control$max_delta) {
      convergence <- FALSE
    }
  } else if (control$delta == control$previous_delta) {
    stop <- TRUE
    reason <- "repeated_delta"
    if (control$delta > control$max_delta) {
      convergence <- FALSE
    }
  } else if (
    control$n_step > 5L &&
      control$delta == source_gllrm_history_value(control, control$n_step - 2L)
  ) {
    stop <- TRUE
    reason <- "two_back_repeated_delta"
    if (control$delta > control$max_delta) {
      convergence <- FALSE
    }
  } else if (control$n_step %% 50L == 0L && control$delta > 10) {
    stop <- TRUE
    reason <- "periodic_50_large_delta"
    if (control$delta > control$max_delta) {
      convergence <- FALSE
    }
  } else if (isTRUE(control$finish) && control$n_finish > 10L) {
    stop <- TRUE
    reason <- "finish_count_plateau"
    convergence <- FALSE
  }

  list(
    stop = stop,
    reason = reason,
    convergence = convergence,
    recurring = recurring
  )
}

#' Record one completed GLLRM IPF discrepancy
#'
#' Mirrors the bookkeeping order at the end of
#' `source/GLLRM_ESTIM.txt::Take_an_IPF_step`, then evaluates the source stop
#' functions.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param control A source GLLRM stopping-control list.
#' @param delta Discrepancy produced by the completed IPF step.
#' @return A list containing updated `control` and its `decision`.
#' @keywords internal
source_gllrm_observe_delta <- function(control, delta) {
  control$n_step <- control$n_step + 1L
  control$finish <- FALSE
  control$previous_delta <- control$delta
  control$delta <- as.numeric(delta)

  if (control$n_step == 1L) {
    control$initial_delta <- control$delta
    control$min_delta <- control$delta
    control$min_delta_step <- 1L
  }
  if (control$delta < control$min_delta) {
    control$min_delta <- control$delta
    control$min_delta_step <- control$n_step
  }
  if (control$previous_delta <= control$delta) {
    control$finish <- TRUE
    control$n_finish <- control$n_finish + 1L
  } else {
    control$finish <- FALSE
    control$n_finish <- 0L
  }
  control$delta_history[as.character(control$n_step)] <- control$delta

  decision <- source_gllrm_stop_decision(control)
  control$convergence <- decision$convergence
  control$recurring <- decision$recurring
  if (isTRUE(decision$stop)) {
    control$stop_reason <- decision$reason
  }
  list(control = control, decision = decision)
}

#' Apply DIGRAM's post-stop convergence acceptance
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param control A stopped source GLLRM control state.
#' @return A single logical value; strict `delta < 0.1` can accept a fit that
#'   stopped with its pre-finalization convergence flag false.
#' @keywords internal
source_gllrm_final_convergence <- function(control) {
  isTRUE(control$convergence) || control$delta < 0.1
}
