#' Summary tables for gRm objects
#'
#' Return R-native summary tables for objects produced by the public gRm
#' API.
#'
#' @name gRm-summary
#' @aliases summary.gRm_analysis summary.gRm_model summary.gRm_fit
#'   summary.gRm_screen
#'   summary.gRm_local_dependence
#'   summary.gRm_dif summary.gRm_global_homogeneity
#'   summary.gRm_m2 summary.gRm_m3 print.summary.gRm
#' @usage
#' \method{summary}{gRm_analysis}(object, ...)
#' \method{summary}{gRm_model}(object, ...)
#' \method{summary}{gRm_fit}(object, which = NULL, ...)
#' \method{summary}{gRm_screen}(object, ...)
#' \method{summary}{gRm_local_dependence}(object, ...)
#' \method{summary}{gRm_dif}(object, ...)
#' \method{summary}{gRm_global_homogeneity}(object, which = NULL, ...)
#' \method{summary}{gRm_m2}(object, ...)
#' \method{summary}{gRm_m3}(object, ...)
#' \method{print}{summary.gRm}(x, ...)
#' @param object An gRm analysis, model, fit, screen, or result object.
#' @param which Optional named summary section to return for summary methods
#'   with multiple views. Available for `gRm_fit` and
#'   `gRm_global_homogeneity`; single-view summaries do not accept `which`.
#' @param x A `summary.gRm` object returned by `summary()`.
#' @param ... Reserved for S3 dispatch compatibility; ignored by all gRm
#'   summary and summary-print methods.
#' @return A `summary.gRm` list. The requested tables are available both by
#'   section name and through the `tables` component. `print.summary.gRm()`
#'   returns the input object invisibly.
#' @details
#' The public output layer is based on `summary()` rather than DIGRAM
#' fixed-width text reports. The returned object contains R data frames suitable
#' for inspection, downstream analysis, or formatting by user code. Only
#' summaries with multiple user-selectable views accept `which`.
#' Printed `gRm_analysis` summaries start with a compact
#' `gRm: Graphical Log-Linear Rasch Model analysis` header showing the data
#' argument label, row count, item and exogenous roles, their constructed
#' project level counts, ID variable, and score groups without repeating the
#' selected data table. It then prints a compact score-group distribution
#' derived from the constructed project score data and analysis cutpoints, with
#' `Score`, `Count`, `Percent`, and `Cumulative` columns. The returned summary
#' object still contains the `"data"` and `"score_groups"` tables for
#' programmatic access. `gRm_analysis` has one summary surface and does not
#' accept `which`.
#'
#' Available sections are class-specific:
#'
#' * `gRm_analysis`: `"data"`, `"score_groups"`.
#' * `gRm_model`: `"model"`. This prints and returns the complete model
#'   specification surface, including the model summary, LD terms, and DIF
#'   terms. Printed model summaries show LD and DIF terms as compact term
#'   lists; the returned summary object retains the full term tables.
#' * `gRm_fit`: `summary(fit)` prints fitted item parameters, thresholds, and
#'   fitted LD/DIF terms. `which = "parameters"` returns and prints only item
#'   parameter summaries; `which = "thresholds"` returns and prints only
#'   threshold rows. Fit status is shown by printing the fitted object.
#' * `gRm_screen`: `summary(screen_obj)` prints one SCREEN J test surface with
#'   final public decision rows for local-dependence tests, DIF tests, and
#'   score effects. These rows are source-backed but are not a complete DIGRAM
#'   oracle report rendering and not a dump of every intermediate SCREEN J
#'   source table. The LD section shows both directed partial-gamma estimates
#'   and p-values, WPG, the directed-gamma sum, and the final model decision.
#'   Local-dependence and DIF candidate evidence use the global
#'   Benjamini-Hochberg FDR 0.05 threshold; stricter global cutoffs for FDR
#'   0.01 and 0.001 are retained in `attr(x, "bh")`. DIGRAM's subsequent
#'   greedy LD procedure may retain negative evidence provisionally, but the
#'   final screen model includes a provisional LD pair only when its two
#'   directed partial gammas have a strictly positive sum. Excluded negative LD
#'   remains visible in the `Decision` column and is never passed to
#'   `gllrm(screen_obj)`. Score-effect rows follow the source score-effect
#'   screening routine, not the global LD/DIF BH table. Final selected model
#'   terms are printed and available as `$selected`, global LD/DIF BH
#'   thresholds as `attr(x, "bh")`, and all final model terms as
#'   `attr(x, "model_terms")`.
#'   The SCREEN J DIF table has separate
#'   `Chisq` / `Pr(>Chisq)` and `Gamma` / `Pr(>|Gamma|)` columns because the
#'   source algorithm chooses the statistic from the exogenous variable shape:
#'   gamma is used when the exogenous variable is binary by category count
#'   (`raw_max == 2`) or ordinal by source variable type (`vtype > 2`), and
#'   chi-square otherwise. The unused statistic-family cells print as blanks
#'   and remain `NA` in the returned table.
#' * `gRm_local_dependence`: `"tests"`. Printed summaries show the complete
#'   local-dependence diagnostic test table directly below the title, with
#'   item names, chi-square statistics, degrees of freedom, `format.pval()`-
#'   formatted p-values, WPG, yes/no convergence status, delta, and a final
#'   blank-headed `*` marker column for rows selected by the
#'   Benjamini-Hochberg threshold. The marker note below the table reports the
#'   FDR value and applied threshold. The returned summary object keeps
#'   selected rows as `$selected` and the BH cutoff metadata as
#'   `attr(x, "bh")`.
#' * `gRm_dif`: `"tests"`. Printed summaries show the complete DIF
#'   candidate test table directly below the title, with item name, exogenous
#'   variable name, chi-square statistic, degrees of freedom,
#'   `format.pval()`-formatted p-value, gamma, yes/no convergence status,
#'   yes/no output-stability status, delta, and a final blank-headed `*`
#'   marker column for rows selected by the Benjamini-Hochberg threshold. The
#'   marker note below the table reports the FDR value and applied threshold.
#'   The returned summary object keeps selected rows as `$selected`, included
#'   DIF rows as `$included`, and the BH cutoff metadata as `attr(x, "bh")`.
#' * `gRm_global_homogeneity`: `summary(result)` prints the source-facing
#'   global homogeneity report: the global test table, fitted score-group
#'   table, item means by score group, and, when the fitted model has LD and/or
#'   DIF terms, uniform local-dependence and uniform DIF interaction sections.
#'   Section-specific calls may use `which = "test"`,
#'   `which = "score_groups"`, `which = "item_means"`,
#'   `which = "uniform_ld"`, or `which = "uniform_dif"`. The uniform sections
#'   show one row per fitted LD/DIF term with observed and expected gamma
#'   columns named by the data-dependent score-group labels, followed by
#'   chi-square, degrees of freedom, and p-value columns. Residual and marker
#'   source-status metadata are available as `attr(x, "source_status")`; item
#'   residual and marker cells are not printed because the recovered Pascal
#'   source does not source-back the exact runtime residual materialization.
#' * `gRm_m2` and `gRm_m3`: one summary surface. Printed summaries show the
#'   source-backed M2/M3 aggregate rows, item-trait interaction aggregate,
#'   invariance decomposition by exogenous variable, and prepared margin rows.
#'   They do not accept `which`; programmatic tables are available from the
#'   returned summary object and from `x$values`.
#'
#' @examples
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
#' model <- gllrm(analysis)
#' summary(analysis)
#' summary(model)
#' @seealso [gRm()], [gllrm()], [fit()], [screen()],
#'   [score_effects()], [item_fit()], [local_dependence()], [dif()],
#'   [global_homogeneity()]
NULL

#' Print gRm objects
#'
#' Print compact status information for public gRm objects.
#'
#' @name gRm-print
#' @aliases print.gRm_analysis print.gRm_model print.gRm_fit
#'   print.gRm_screen
#'   print.gRm_direct_table
#'   print.gRm_local_dependence print.gRm_dif
#'   print.gRm_global_homogeneity print.gRm_m2 print.gRm_m3
#' @usage
#' \method{print}{gRm_analysis}(x, ...)
#' \method{print}{gRm_model}(x, ...)
#' \method{print}{gRm_fit}(x, ...)
#' \method{print}{gRm_screen}(x, ...)
#' \method{print}{gRm_direct_table}(x, ...)
#' \method{print}{gRm_local_dependence}(x, ...)
#' \method{print}{gRm_dif}(x, ...)
#' \method{print}{gRm_global_homogeneity}(x, ...)
#' \method{print}{gRm_m2}(x, ...)
#' \method{print}{gRm_m3}(x, ...)
#' @param x An gRm analysis, model, fit, screen, or result object.
#' @param ... Reserved for S3 dispatch compatibility; ignored by all gRm
#'   print methods.
#' @return The input object, invisibly.
#' @details
#' The print methods are intentionally compact. They identify the object type
#' and key counts or status fields, but they do not expose the full numerical
#' result surface. Analysis print output shows the source label, row count, ID,
#' variable counts, and score-group count. Model print output intentionally
#' matches the complete human-facing model summary text because a model
#' specification is compact. Fitted-model print output shows convergence,
#' iteration, delta, and log-likelihood status. SCREEN J print output shows
#' compact model-search status, including inference mode, source-backed tested
#' relation count, selected LD/DIF term counts, and selected score-effect
#' count. Local-dependence and DIF diagnostic print methods show candidate
#' counts, Benjamini-Hochberg selected counts, convergence status counts, and
#' the applied BH threshold; DIF print output also reports included DIF terms
#' and output-stable fits. Direct table
#' results, such as those returned by `item_fit()` and `score_effects()`, print
#' their selected table directly with a `gRm:` title. Use `summary()` for
#' R-native summary result tables where those summary methods are available.
#' @examples
#' data <- data.frame(
#'   ID = 1:4,
#'   I1 = c(0, 1, 0, 1),
#'   I2 = c(1, 0, 1, 0)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), id = "ID", score_cuts = c(1L, 2L))
#' print(analysis)
#' print(gllrm(analysis))
#' @seealso [summary()], [gRm()], [gllrm()], [fit()]
NULL
