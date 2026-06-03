#' Run source-faithful DIGRAM SCREEN J
#'
#' Run DIGRAM SCREEN J model discovery from an analysis object.
#'
#' @usage screen(project, inference = c("asymptotic", "exact", "repeated"),
#'   nsim = 1000L, seed = 9L, critlevel = NULL, risk = NULL, jobs = 1L, ...)
#' @param project A `gRm_analysis` object from [gRm()] or
#'   [read_digram_project()], or a source-shaped `gRm_project`.
#' @param inference Inference mode. `"asymptotic"` maps to DIGRAM command 17,
#'   `"exact"` to command 2, and `"repeated"` to command 74.
#' @param nsim Number of simulations for exact branches.
#' @param seed Random seed for exact branches.
#' @param critlevel Optional DIGRAM exact command critical level on the
#'   per-1000 scale.
#' @param risk Optional DIGRAM repeated-Monte-Carlo risk on the per-1000 scale.
#' @param jobs Reserved for future parallel screen implementations.
#' @param ... Reserved for S3 dispatch compatibility; no additional SCREEN J
#'   arguments are currently supported.
#' @return A `gRm_screen` object containing SCREEN J test values, selected
#'   model terms, score-effect rows, Benjamini-Hochberg thresholds, and
#'   exact-command metadata when relevant.
#' @details
#' SCREEN J evaluates candidate LD, DIF, and score-effect relations using the
#' source-shaped SCREEN J algorithm. The result can be inspected with
#' `summary(screen_obj, which = "tests")`, `"selected"`, `"all"`,
#' `"score_effects"`, or `"bh"`. Passing a screen object to [gllrm()] creates
#' a model from selected LD and DIF terms.
#'
#' The `inference` modes map to DIGRAM command-state conventions:
#' `"asymptotic"` uses the asymptotic path, `"exact"` uses fixed Monte Carlo
#' exact inference, and `"repeated"` uses the repeated Monte Carlo command.
#' `nsim`, `seed`, `critlevel`, and `risk` configure the exact command state.
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   group = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "group", id = "ID")
#' scr <- screen(analysis, inference = "asymptotic")
#' summary(scr, which = "selected")
#' model <- gllrm(scr)
#' }
#' @seealso [gRm()], [gllrm()], [score_effects()]
#' @export
screen <- function(project,
                   ...) {
  UseMethod("screen")
}

#' @export
screen.gRm_analysis <- function(project,
                   inference = c("asymptotic", "exact", "repeated"),
                   nsim = 1000L,
                   seed = 9L,
                   critlevel = NULL,
                   risk = NULL,
                   jobs = 1L,
                   ...) {
  analysis <- project
  inference <- match.arg(inference)
  exact_state <- gRm_exact_command_state_public(
    inference,
    nsim = nsim,
    seed = seed,
    critlevel = critlevel,
    risk = risk
  )
  values <- screen_j_values(
    analysis$project,
    exact = exact_state$exact,
    repeated = exact_state$sequential,
    nsim = exact_state$nsim,
    seed = exact_state$seed,
    exact_state = exact_state
  )
  out <- list(
    analysis = analysis,
    project = analysis$project,
    values = values,
    inference = inference,
    exact_state = exact_state,
    nsim = exact_state$nsim,
    seed = if (exact_state$exact) exact_state$seed else NA_integer_,
    jobs = as.integer(jobs %||% 1L),
    terms = NULL,
    source_trace = c(values$source_status %||% character(), api = "screen"),
    unmodeled = character(),
    warnings = character(),
    call = match.call()
  )
  class(out) <- c("gRm_screen", "list")
  out$terms <- model_terms(out)
  out
}

#' @export
screen.gRm_item_analysis <- screen.gRm_analysis

#' @export
screen.gRm_project <- function(project,
                                  inference = c("asymptotic", "exact", "repeated"),
                                  nsim = 1000L,
                                  seed = 9L,
                                  critlevel = NULL,
                                  risk = NULL,
                                  jobs = 1L,
                                  ...) {
  screen.gRm_analysis(
    as_gRm_analysis(project),
    inference = inference,
    nsim = nsim,
    seed = seed,
    critlevel = critlevel,
    risk = risk,
    jobs = jobs
  )
}

screen.default <- function(project, ...) {
  stop("screen() requires an analysis returned by gRm or read_digram_project.", call. = FALSE)
}

#' @export
print.gRm_screen <- function(x, ...) {
  cat("<gRm_screen>\n")
  cat("  inference: ", x$inference, "\n", sep = "")
  cat("  proposed LD terms: ", nrow(x$terms$ld), "\n", sep = "")
  cat("  proposed DIF terms: ", nrow(x$terms$dif), "\n", sep = "")
  invisible(x)
}
