#' Run source-faithful DIGRAM SCREEN J
#'
#' Run DIGRAM SCREEN J model discovery from an analysis object.
#'
#' @usage
#' screen(project, ...)
#' \method{screen}{gRm_analysis}(project,
#'   inference = c("asymptotic", "exact", "repeated"), nsim = 1000L,
#'   seed = 9L, critlevel = NULL, risk = NULL, ...)
#' @aliases screen.gRm_analysis
#' @param project A `gRm_analysis` object from [gRm()] or
#'   [read_digram_project()], or a source-shaped `gRm_project`.
#' @param inference Inference mode. `"asymptotic"` maps to DIGRAM command 17,
#'   `"exact"` to command 2, and `"repeated"` to command 74.
#' @param nsim Single non-negative integer-like number of simulations for exact
#'   branches. The source convention `0` requests the DIGRAM default of `1000`.
#' @param seed Single integer-like random seed for exact branches.
#' @param critlevel Optional DIGRAM exact command critical level on the
#'   per-1000 scale; when supplied it must be a single non-negative
#'   integer-like value.
#' @param risk Optional DIGRAM repeated-Monte-Carlo risk on the per-1000 scale;
#'   when supplied it must be a single non-negative integer-like value.
#' @param ... Reserved for S3 dispatch compatibility; no additional SCREEN J
#'   arguments are currently supported.
#' @return A `gRm_screen` object containing SCREEN J test values, selected
#'   model terms, score-effect rows, Benjamini-Hochberg thresholds, and
#'   exact-command metadata when relevant.
#' @details
#' SCREEN J evaluates candidate LD, DIF, and score-effect relations using the
#' source-shaped SCREEN J algorithm. `print(screen_obj)` shows compact status
#' information; `summary(screen_obj)` prints the final public SCREEN J decision
#' rows for local dependence, DIF, and score effects. The summary is a
#' human-facing test surface, not a complete DIGRAM oracle report rendering and
#' not a dump of all intermediate source calculations. It deliberately omits
#' internal partial item-pair matrices, initial item-bias matrices,
#' spurious-DIF reduction iterations, multiple-DIF reduction iterations, and
#' score-effect intermediate screening iterations from the printed output.
#' Those source-shaped values remain in `screen_obj$values` for validation,
#' debugging, and source-faithfulness audits.
#'
#' Printed SCREEN J summaries show both directed partial-gamma estimates and
#' p-values for each unordered item pair, together with WPG, the sum of the two
#' directed gammas, and the final decision. Local-dependence and DIF candidate
#' evidence use the global Benjamini-Hochberg FDR 0.05 threshold; stricter
#' global cutoffs for FDR 0.01 and 0.001 are retained in
#' `attr(summary(screen_obj), "bh")`. The source then applies a greedy
#' local-dependence evidence procedure. Its provisional negative-LD rows remain
#' visible with decision `"negative LD; not included"`, but—following DIGRAM's
#' version 3.37 rule—only a provisional pair whose two directed partial gammas
#' have a strictly positive sum enters the screen model. Thus `*` marks a final
#' included LD term rather than merely one directed p-value below the BH
#' threshold. Score-effect rows follow the source score-effect screening
#' routine, not the global LD/DIF BH table. The returned summary prints and
#' retains the final model terms as `$selected`, keeps global LD/DIF BH
#' thresholds as `attr(summary(screen_obj), "bh")`, and retains all final model
#' terms as `attr(summary(screen_obj), "model_terms")`.
#' In the printed DIF table, `Chisq` / `Pr(>Chisq)` and `Gamma` /
#' `Pr(>|Gamma|)` are separate statistic families. SCREEN J uses the gamma
#' statistic when the exogenous variable is binary by category count
#' (`raw_max == 2`) or ordinal by source variable type (`vtype > 2`) in the
#' source-shaped background metadata. Otherwise it uses the chi-square
#' statistic. Cells for the statistic family not used for a row are printed
#' blank; the returned summary table keeps those cells as `NA`.
#'
#' Passing a screen object to [gllrm()] creates a model from final selected LD
#' and DIF terms; provisional negative LD is never passed to the model.
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
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )
#' scr <- screen(analysis, inference = "asymptotic")
#' summary(scr)
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
                   ...) {
  reject_public_dots(...)
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
screen.gRm_project <- function(project,
                                  inference = c("asymptotic", "exact", "repeated"),
                                  nsim = 1000L,
                                  seed = 9L,
                                  critlevel = NULL,
                                  risk = NULL,
                                  ...) {
  reject_public_dots(...)
  screen.gRm_analysis(
    as_gRm_analysis(project),
    inference = inference,
    nsim = nsim,
    seed = seed,
    critlevel = critlevel,
    risk = risk
  )
}

#' @export
screen.default <- function(project, ...) {
  stop("screen() requires an analysis returned by gRm or read_digram_project.", call. = FALSE)
}

#' @export
print.gRm_screen <- function(x, ...) {
  tables <- public_screen_summary_tables(x)
  bh <- x$values$bh %||% list()
  exact_state <- x$exact_state %||% list()
  cat("gRm: SCREEN J tests\n\n")
  cat("  Inference: ", x$inference %||% "asymptotic", "\n", sep = "")
  if (isTRUE(exact_state$exact)) {
    cat("  Simulations: ", summary_scalar(x$nsim %||% exact_state$nsim %||% NA_integer_), "\n", sep = "")
    if (isTRUE(exact_state$sequential)) {
      cat("  Sequential limit: ", summary_scalar(exact_state$seq_limit %||% NA_integer_), "\n", sep = "")
    }
    cat("  Seed: ", summary_scalar(x$seed %||% exact_state$seed %||% NA_integer_), "\n", sep = "")
  }
  cat("  Tested relations: ", summary_scalar(bh$n_tests %||% NA_integer_), "\n", sep = "")
  cat("  Selected LD terms: ", summary_scalar(nrow(tables$selected_ld)), "\n", sep = "")
  cat("  Selected DIF terms: ", summary_scalar(nrow(tables$selected_dif)), "\n", sep = "")
  cat("  Selected score effects: ", summary_scalar(sum(tables$score_effects[[" "]] %in% "*")), "\n\n", sep = "")
  cat("Use summary(x) to show the SCREEN J test tables.\n")
  invisible(x)
}
