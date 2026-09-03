#' CM2 fit diagnostic
#'
#' Compute DIGRAM's source-backed CM2 diagnostic for a fitted gRm model.
#'
#' @param fit A fitted `gRm_fit` object returned by [fit()]. An unfitted
#'   `gRm_analysis`, a `gRm_model`, raw data, and other object types are not
#'   accepted. Both fitted Rasch and fitted GLLRM specifications are supported.
#'   CM2/CM3 use the supplied fit's item and any LD/DIF parameters together
#'   with its complete fitted exogenous-variable scope; the fit is not modified.
#' @param items Focal items whose diagnostic margins are to be tested. The
#'   permitted values are:
#'
#'   * `NULL` (the default), which selects every fitted item and corresponds to
#'     leaving DIGRAM's item-selection prompt blank;
#'   * a character vector containing at least two distinct, exact,
#'     case-sensitive R-facing fitted item names; or
#'   * a numeric or integer vector containing at least two distinct, finite,
#'     whole, one-based fitted-item indices.
#'
#'   Character and numeric selectors cannot be mixed. The function rejects
#'   zero-length vectors, duplicates, missing or non-finite values, fractional
#'   or out-of-range indices, empty names, factors, logical vectors, matrices,
#'   arrays, exogenous-variable names, and the score-group marker. The literal
#'   string `"all"` has no special meaning. Input order does not control output
#'   order: valid selections are normalized to fitted/source item order because
#'   DIGRAM stores the selection as a Boolean set. Explicitly naming or indexing
#'   every item is numerically equivalent to `NULL`, but the recorded selection
#'   mode remains explicit. Selection affects diagnostic margins only, as
#'   described under **Selection and score-group scope**.
#' @param score_cuts Total-score upper bounds used to construct the temporary
#'   score-group variable for this diagnostic. `NULL` (the default) uses the
#'   resolved cuts stored in the fitted analysis; if fewer than two stored cuts
#'   are available, the source-faithful automatic cuts are recomputed. An
#'   explicit value must be a numeric or integer vector of at least two finite,
#'   non-missing whole values, in strictly increasing order, within
#'   `0..sum(maximum item scores)`. After clipping to DIGRAM's observed,
#'   non-extreme score range, the cuts must define at least two usable groups.
#'   Each value is an inclusive upper boundary and subsequent groups begin one
#'   score above the preceding boundary. The returned `score_groups` bounds are
#'   display-clipped to DIGRAM's observed, non-extreme range. For the actual
#'   CM2/CM3 table counts, the source `Scoregruppe` rule extends the first
#'   retained group down to score zero and the last retained group through the
#'   highest possible total score, so both extreme scores remain classified.
#'   Unlike the `score_cuts` argument to [gRm()], the string `"auto"` is not
#'   accepted here. Explicit cuts are local to the returned diagnostic: they
#'   neither mutate the analysis nor refit the model.
#' @param bootstrap A single non-missing logical value. `FALSE` (the default)
#'   returns only the deterministic/asymptotic diagnostic and corresponds to
#'   zero simulations at DIGRAM 7.04's prompt. `TRUE` additionally requests the
#'   source-shaped parametric bootstrap. If the fitted model exceeds a preserved
#'   DIGRAM generator bound, the ordinary diagnostic is still returned and
#'   `result$values$bootstrap$possible` is `FALSE`; no approximate fallback is
#'   used and no bootstrap p-values are fabricated.
#' @param nsim A single finite, whole numeric or integer value from `1` through
#'   `.Machine$integer.max`: the number of bootstrap samples to generate when
#'   `bootstrap = TRUE`. The default is `100L`. This is the attempted sample
#'   count; `result$values$bootstrap$nused` is the possibly smaller number whose
#'   refits satisfy DIGRAM's strict acceptance rule. The value is validated but
#'   not used when `bootstrap = FALSE`.
#' @param seed Bootstrap random-stream initialization. Use `NULL` (the default)
#'   to reproduce Delphi 4 `Randomize`, which initializes from UTC milliseconds
#'   since midnight, or supply one finite whole numeric/integer value in
#'   `0..4294967295` as the unsigned 32-bit Delphi `RandSeed` bit pattern. An
#'   explicit value gives an exactly reproducible stream; values above
#'   `2147483647` correspond to the negative `LongInt` values that Delphi may
#'   display but must be supplied to R in unsigned form. The bootstrap uses a
#'   private RNG and never reads or changes R's global `.Random.seed`. The value
#'   is validated but not used when `bootstrap = FALSE`.
#' @param reestimate A single non-missing logical value controlling each
#'   generated sample. `TRUE` (the default) starts a fresh fit of the same full
#'   GLLRM in every sample and is the canonical DIGRAM 7.04 command-197
#'   behavior. `FALSE` evaluates every generated sample under the supplied
#'   fitted parameters and is an explicit R extension. This argument concerns
#'   bootstrap-sample refitting only; it never refits or mutates `fit` itself.
#'   In the `FALSE` branch, every sample reuses the observed fit's source
#'   stopping-step `report_delta`; consequently all generated samples pass the
#'   strict `< 0.1` acceptance gate when that value is below `0.1`, or none do
#'   otherwise.
#'   DIGRAM's final "Item parameters were reestimated" footer is unconditional
#'   and therefore does not identify this branch. The value is validated but
#'   not used when `bootstrap = FALSE`.
#' @param bootstrap_max_step A single finite, whole numeric or integer value
#'   from `1` through `.Machine$integer.max`, giving the maximum fitting steps for each
#'   bootstrap refit. The default, `5000L`, is hard-coded in
#'   `SKbias8.Estimate_the_GLLRM`; other values are an explicit R runtime-control
#'   extension. DIGRAM's separate `Nsteps = 500` generator argument is unused
#'   by the preserved source. This control has an operational effect only when
#'   both `bootstrap = TRUE` and `reestimate = TRUE`, but it is validated on
#'   every call.
#' @param bootstrap_jobs `NULL`, or a single finite, whole numeric/integer value
#'   from `1` through `.Machine$integer.max`, giving the requested maximum number of
#'   workers for bootstrap refitting and diagnostic analysis. `NULL` and the
#'   default `1L` both request serial execution. On POSIX systems, larger values
#'   use ordered fork workers and are capped by `nsim` and the detected logical
#'   cores; Windows uses one worker. Sample generation always remains serial so
#'   worker scheduling cannot change the Delphi RNG stream or replicate order.
#'   The value is validated but not used when `bootstrap = FALSE`.
#' @param keep_bootstrap_samples A single non-missing logical value. When
#'   `bootstrap = TRUE`, `TRUE` retains every generated zero-based item-response
#'   matrix and one-based exogenous matrix in
#'   `result$values$bootstrap$samples`; `FALSE` (the default) returns an empty
#'   `samples` list to reduce returned-object size. Whenever a supported
#'   simulation runs, per-replicate fit, aggregate, and margin statistic tables
#'   are retained under either setting. Retaining all response matrices can use
#'   substantial memory. Sample generation still requires working memory, so
#'   `FALSE` controls result retention rather than eliminating generation-time
#'   memory use. The value is validated but not used when `bootstrap = FALSE`.
#' @param resample_score_distribution A single non-missing logical value for
#'   Pascal's `ResampleScoreDist` control. The preserved source branch is empty,
#'   so both `FALSE` (the default) and `TRUE` keep the observed score
#'   distribution within complete exogenous strata. `TRUE` does **not** request
#'   a different resampling algorithm; the requested value and the source no-op
#'   mode are retained in bootstrap metadata. The value is validated but not
#'   used when `bootstrap = FALSE`.
#' @param ... Reserved for future extensions. It must be empty; every named or
#'   unnamed additional argument is rejected.
#' @return A `gRm_cm2` list with the standard gRm result components
#'   (`analysis`, `project`, `fit`, `metadata`, `call`, and related provenance)
#'   and the following CM2 results under `result$values`:
#'
#'   * `tests`: one row per prepared two-way margin, in DIGRAM source order.
#'     `margin_type` is one of `"item_item"`, `"item_exogenous"`, or
#'     `"item_score_group"`; the row contains public variable names, source
#'     labels, `chi_square`, `degrees_of_freedom`, and the asymptotic `p_value`.
#'   * `aggregate`: the CM2 sum, with its summed degrees of freedom and source
#'     chi-square-tail p-value; `item_trait` and `invariance` contain the
#'     corresponding source decompositions.
#'   * `cm2_bh`: the 5%, 1%, and 0.1% observed-data Benjamini--Hochberg critical
#'     p-values and the source 85-p-value capacity metadata.
#'   * `selected_items` and `selection`: the resolved source-order item table
#'     and the complete selection record. `selection$mode` is
#'     `"default_all"`, `"explicit_names"`, or `"explicit_indices"`;
#'     `requested_items` preserves the evaluated request, while
#'     `resolved_items`, `exogenous_scope`, and `score_total_scope` describe the
#'     actual diagnostic scope. `metadata$selection` is the same record, and
#'     `metadata$items` and `metadata$item_indices` are convenience vectors.
#'   * `score_cuts` and `score_groups`: the resolved cuts and retained inclusive
#'     group bounds; `n_two_way_margins` records the prepared row count.
#'   * `bootstrap`: bootstrap status, controls, capability, and results. See
#'     **Bootstrap result states** for the three possible structures.
#'
#'   When a supported bootstrap is run, `tests`, `aggregate`, `item_trait`, and
#'   `invariance` also receive `bootstrap_extreme_count`, `bootstrap_nused`, and
#'   `bootstrap_p_value` columns. Use `summary(result)` for presentation tables;
#'   use `result$values` for programmatic access.
#' @details
#' `cm2()` is a diagnostic of the current fitted model. It does not search for,
#' select, add, or refit local-dependence or DIF terms.
#'
#' In the DIGRAM source, CM2 is not a separately registered command. It is the
#' two-way-margin part of DIGRAM's CM3 command. The source `Prepare_CM3tests`
#' routine prepares item-item margins, item-exogenous margins, and
#' item-score-group margins for the selected items. Item-item margins that are
#' already included as local-dependence terms are skipped. DIGRAM 7.04 corrected
#' the item-exogenous filter to `ItemBias(.i1,i2.)`, so modeled DIF margins are
#' skipped in the documented item-first, exogenous-second coordinate order.
#' Item-score-group rows are included for selected items.
#'
#' The CM2 aggregate is the sum of the prepared two-way Pearson chi-square
#' contributions. Degrees of freedom and p-values follow the DIGRAM source,
#' using the package's `source_pfchi()` implementation for chi-square tail
#' probabilities. Public output uses R variable names. DIGRAM's compact
#' alphabetic labels are retained internally for source tracing and oracle
#' validation.
#'
#' @section Selection and score-group scope:
#' `items` controls only the focal diagnostic margins. Every fitted exogenous
#' variable remains automatic, and score groups are formed from the total over
#' all fitted items, including unselected items. Fitted LD and DIF terms,
#' including LD partners that are not selected, remain in current-model
#' probabilities. Active terms suppress only their eligible CM2 item-item or
#' item-exogenous rows. The observed fit is never mutated or refitted.
#'
#' `score_cuts` changes only the temporary score-group rows in this result.
#' Score-group membership uses the zero-based total over all fitted items, not
#' the subtotal over `items`. The returned display bounds are clipped to
#' DIGRAM's observed non-extreme score window, whereas the source table-counting
#' lookup extends the first and last retained groups to scores zero and the
#' highest possible score. Complete records at both extremes therefore remain
#' eligible for CM2/CM3 score-group margins.
#'
#' Bootstrap generation and optional refitting also use the complete fitted
#' model. The selected items and diagnostic score cuts control only the fixed
#' margin list evaluated in the observed data and in every accepted sample.
#'
#' The R diagnostics use dynamically sized margin lists and do not impose
#' DIGRAM command 197's static `DATALENGTH = 1600` combined CM2/CM3 storage
#' bound. A result whose prepared margin list exceeds that executable limit is
#' a supported R computation, but it lies outside direct DIGRAM
#' executable-oracle capacity and is not claimed as executable-oracle
#' validated.
#'
#' @section Bootstrap algorithm and calibration:
#' With `bootstrap = TRUE`, score counts are fixed within complete exogenous
#' strata and response patterns are generated conditionally from the supplied
#' fitted GLLRM. Active LD and DIF parameters participate in the same component
#' probabilities as fitting. Generation uses Delphi 4's wrapped 32-bit linear
#' congruential `Random` stream, preserves both occurrences of source
#' `x := 1 - random`, and follows DIGRAM's component-score and within-component
#' response-pattern draw order. When `reestimate = TRUE`, refits start from
#' fresh source parameter values and use `bootstrap_max_step`; the acceptance
#' discrepancy is Pascal's global `delta` at loop exit, retained as R
#' `report_delta`, rather than the R fitter's later recomputed sufficient-count
#' discrepancy. When `reestimate = FALSE`, the supplied observed fit's
#' `report_delta` is reused for every sample. Source refits use the
#' inclusive score window from zero through the highest possible score, so
#' both deterministic extreme-score blocks enter every fitted margin. The
#' returned `reestimate` metadata always describes the per-bootstrap-sample
#' `CM3_analysis.DoSomething` branch. It is deliberately not inferred from
#' DIGRAM's unconditional cleanup footer in `Finish_random_Gllrm`.
#'
#' A sample contributes only when its source acceptance discrepancy is strictly
#' less than `0.1`. Bootstrap p-values are
#' `sum(p_simulated <= p_observed) / nused` over those accepted samples, with no
#' plus-one correction. Therefore `nsim` is the requested/generated count and
#' `nused` is the calibration denominator. If no sample is accepted,
#' `nused = 0` and every bootstrap p-value is `NA`. `bootstrap_jobs` changes
#' only the scheduling of the RNG-free refit/analysis phase; requested and
#' effective worker counts and the execution mode are recorded in the result.
#'
#' @section Bootstrap capability limits:
#' Generation is available only inside the canonical fixed-array bounds: at
#' least one and at most 15,000 complete item/exogenous records; at most 10
#' fitted items; at most four items in an LD component; at most four multi-item
#' components; and item scores no greater than seven. The Cartesian product of
#' the category counts of all fitted exogenous variables may not exceed 216;
#' the corresponding product for exogenous variables in active DIF terms may
#' not exceed 64. These are full possible-stratum counts, including category
#' combinations absent from the data. A fitted component score may have at most
#' 255 positive response patterns because Pascal stores that count in a `Byte`.
#' An out-of-bound request returns explicit capability reasons without running
#' a distributionally different fallback.
#'
#' @section Bootstrap result states:
#' `result$values$bootstrap` always records one of three states:
#'
#' * Not requested: `enabled = FALSE`, `possible = NA`, `nsim = 0`, and
#'   `nused = 0`.
#' * Requested but outside a source bound: `enabled = TRUE`,
#'   `possible = FALSE`, `source_status = "bootstrap_not_possible"`, and
#'   `reasons` plus `capability$bounds` and `capability$observed` explain the
#'   refusal. The ordinary/asymptotic CM2 or CM3 result remains valid and is
#'   returned normally.
#' * Completed: `enabled = TRUE`, `possible = TRUE`, and the component records
#'   the effective seed and RNG state, requested and accepted sample counts,
#'   fitting and worker controls, source capability, score-by-exogenous groups,
#'   per-replicate fit/aggregate/margin audit tables, aggregate and margin
#'   summaries, and optional generated samples. The calibrated p-values are
#'   also attached to the corresponding public diagnostic tables.
#'
#' Every bootstrap-related argument is type/range validated even when
#' `bootstrap = FALSE`. In that case it has no computational effect and the
#' not-requested state intentionally records `nsim = 0` rather than the unused
#' `nsim` argument.
#'
#' @section Multiple-testing and record handling:
#' DIGRAM 7.04 also applies separate Benjamini-Hochberg evaluations to the
#' observed two-way and three-way margin p-values. Their 5%, 1%, and 0.1% FDR
#' critical levels are retained in `result$values$cm2_bh` and, for [cm3()],
#' `result$values$cm3_bh`. As in the Pascal `Rvector` path, at most the first 85
#' p-values in each block participate.
#'
#' CM2/CM3 record eligibility is independent of the CML estimation score
#' window. Complete score-zero records contribute to observed and expected
#' diagnostic margins. The historical Pascal `Count_IJK` path also retains an
#' item-complete record when exogenous data are missing, while the corresponding
#' expected table and every other margin require complete exogenous values.
#'
#' `summary(result)$tables` is the stable presentation interface. It contains
#' `aggregates`, `item_trait`, `invariance`, `margins`, and `bh`; when bootstrap
#' was requested it also contains `bootstrap`. The `margins` table remains in
#' DIGRAM source order. Completed bootstrap results add the calibrated p-value
#' and accepted-sample count to the applicable presentation tables.
#' @seealso [cm3()], [gllrm()], [fit()], [summary()]
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(1, 2, 1, 2, 1, 2, 2, 1),
#'   I2 = c(2, 1, 2, 1, 2, 1, 2, 1),
#'   I3 = c(1, 1, 2, 1, 2, 2, 1, 2)
#' )
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), id = "ID")
#' fit0 <- fit(gllrm(analysis))
#' cm2_result <- cm2(
#'   fit0,
#'   items = c("I1", "I3"),
#'   score_cuts = c(1L, 3L)
#' )
#' cm2_by_index <- cm2(fit0, items = c(1L, 3L))
#' summary(cm2_result)
#' cm2_result$metadata$selection
#' head(cm2_result$values$tests)
cm2 <- function(fit,
               items = NULL,
               score_cuts = NULL,
               bootstrap = FALSE,
               nsim = 100L,
               seed = NULL,
               reestimate = TRUE,
               bootstrap_max_step = 5000L,
               bootstrap_jobs = 1L,
               keep_bootstrap_samples = FALSE,
               resample_score_distribution = FALSE,
               ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  bootstrap_control <- normalize_cm2_cm3_bootstrap_control(
    bootstrap,
    nsim,
    seed,
    reestimate,
    bootstrap_max_step,
    bootstrap_jobs,
    keep_bootstrap_samples,
    resample_score_distribution
  )
  values <- cm2_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    bootstrap_control = bootstrap_control
  )
  new_cm2_cm3_result(fit, values, "gRm_cm2", "cm2", match.call())
}

#' CM3 fit diagnostic
#'
#' Compute DIGRAM's source-backed CM3 diagnostic for a fitted gRm model.
#'
#' @inheritParams cm2
#' @return A `gRm_cm3` list with the standard gRm result components
#'   (`analysis`, `project`, `fit`, `metadata`, `call`, and related provenance)
#'   and the following CM3 results under `result$values`:
#'
#'   * `tests`: every prepared two-way and three-way margin, in DIGRAM source
#'     order. `is_cm2` distinguishes the two-way block. In addition to the CM2
#'     types documented in [cm2()], `margin_type` can be
#'     `"item_item_item"`, `"item_item_exogenous"`,
#'     `"item_item_score_group"`, `"item_exogenous_exogenous"`, or
#'     `"item_exogenous_score_group"`. Rows contain public variable names,
#'     source labels, `chi_square`, `degrees_of_freedom`, and asymptotic
#'     `p_value`.
#'   * `cm2`: the aggregate over the two-way block; `cm3`: the cumulative
#'     aggregate over the complete CM2 block plus every three-way row. These
#'     are not independent or disjoint aggregate tests. `item_trait` and
#'     `invariance` contain the source decompositions.
#'   * `cm2_bh` and `cm3_bh`: separate 5%, 1%, and 0.1% observed-data
#'     Benjamini--Hochberg critical p-values for the two-way and three-way
#'     blocks, including each block's source 85-p-value capacity metadata.
#'   * `selected_items`, `selection`, `score_cuts`, `score_groups`,
#'     `n_two_way_margins`, and `n_three_way_margins`: the resolved diagnostic
#'     scope and its row counts. `metadata$selection`, `metadata$items`, and
#'     `metadata$item_indices` provide the same selection information in a
#'     convenient top-level metadata form.
#'   * `bootstrap`: bootstrap status, controls, capability, and results, with
#'     the same three states described under **Bootstrap result states**.
#'
#'   When a supported bootstrap is run, `tests`, `cm2`, `cm3`, `item_trait`,
#'   and `invariance` also receive `bootstrap_extreme_count`,
#'   `bootstrap_nused`, and `bootstrap_p_value` columns. Use `summary(result)`
#'   for stable presentation tables and `result$values` for the full numerical
#'   payload.
#' @details
#' `cm3()` is a diagnostic for the current fitted model. It does not search for,
#' select, add, or refit local-dependence or DIF terms.
#'
#' `cm3()` implements the deterministic/asymptotic output of DIGRAM's CM3
#' command for the selected items. DIGRAM's CM3 output includes CM2 first; the
#' CM3 aggregate is therefore not a three-way-only statistic. It includes the
#' complete CM2 two-way aggregate and then adds the prepared three-way margin
#' contributions.
#'
#' The three-way rows are prepared after the CM2 rows in Pascal source order:
#' item-item-item rows, item-item-exogenous and item-item-score-group rows, and
#' then item-exogenous-exogenous and item-exogenous-score-group rows. Included
#' two-way local-dependence or DIF terms do not remove related three-way rows;
#' this is source behavior from `Prepare_CM3tests`.
#'
#' Item selection, public variable names, score-group handling, and bootstrap
#' scope are the same as for [cm2()]. Public output uses R variable names, while
#' DIGRAM alphabetic labels remain internal for source/oracle matching.
#'
#' @inheritSection cm2 Selection and score-group scope
#' @inheritSection cm2 Bootstrap algorithm and calibration
#' @inheritSection cm2 Bootstrap capability limits
#' @inheritSection cm2 Bootstrap result states
#' @inheritSection cm2 Multiple-testing and record handling
#' @seealso [cm2()], [gllrm()], [fit()], [summary()]
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(1, 2, 1, 2, 1, 2, 2, 1),
#'   I2 = c(2, 1, 2, 1, 2, 1, 2, 1),
#'   I3 = c(1, 1, 2, 1, 2, 2, 1, 2)
#' )
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), id = "ID")
#' fit0 <- fit(gllrm(analysis))
#' cm3_result <- cm3(fit0, items = c(1L, 3L))
#' cm3_by_name <- cm3(fit0, items = c("I1", "I3"))
#' summary(cm3_result)
#' if (interactive()) {
#'   cm3_bootstrap <- cm3(
#'     fit0,
#'     items = c("I1", "I3"),
#'     bootstrap = TRUE,
#'     nsim = 1000L,
#'     seed = 9
#'   )
#' }
cm3 <- function(fit,
               items = NULL,
               score_cuts = NULL,
               bootstrap = FALSE,
               nsim = 100L,
               seed = NULL,
               reestimate = TRUE,
               bootstrap_max_step = 5000L,
               bootstrap_jobs = 1L,
               keep_bootstrap_samples = FALSE,
               resample_score_distribution = FALSE,
               ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  bootstrap_control <- normalize_cm2_cm3_bootstrap_control(
    bootstrap,
    nsim,
    seed,
    reestimate,
    bootstrap_max_step,
    bootstrap_jobs,
    keep_bootstrap_samples,
    resample_score_distribution
  )
  values <- cm3_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    bootstrap_control = bootstrap_control
  )
  new_cm2_cm3_result(fit, values, "gRm_cm3", "cm3", match.call())
}

#' Internal new cm2 cm3 result helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param fit Fitted gRm model.
#' @param values Values to validate or transform.
#' @param class Internal `class` value used by this helper.
#' @param result Result value to assemble or transform.
#' @param call Captured R call.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_cm2_cm3_result <- function(fit, values, class, result, call) {
  analysis <- fit$analysis %||% fit$spec$analysis
  selected <- values$selected_items %||% data.frame()
  metadata <- list(
    items = selected$item_name %||% character(),
    item_indices = selected$item_index %||% integer(),
    selection = values$selection,
    score_cuts = values$score_cuts %||% integer(),
    bootstrap = values$bootstrap[c(
      "enabled", "possible", "nsim", "nused", "seed", "reestimate",
      "requested_bootstrap_jobs", "bootstrap_jobs", "execution_mode",
      "acceptance_delta"
    )]
  )
  new_gRm_result(
    class = class,
    analysis = analysis,
    fit = fit,
    values = values,
    result = result,
    metadata = metadata,
    call = call
  )
}

#' @export
print.gRm_cm2 <- function(x, ...) {
  reject_public_dots(...)
  print_cm2_cm3_result(x, include_cm3 = FALSE)
}

#' @export
print.gRm_cm3 <- function(x, ...) {
  reject_public_dots(...)
  print_cm2_cm3_result(x, include_cm3 = TRUE)
}

#' Internal print cm2 cm3 result helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param include_cm3 Internal `include_cm3` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_cm2_cm3_result <- function(x, include_cm3) {
  values <- x$values %||% list()
  title <- if (isTRUE(include_cm3)) "gRm: CM3 fit diagnostic" else "gRm: CM2 fit diagnostic"
  cat(title, "\n\n", sep = "")
  selection <- cm2_cm3_selection_metadata(values)
  cat("  Selection mode: ", selection$mode, "\n", sep = "")
  cat("  Selection: ", cm2_cm3_selection_label(values), "\n", sep = "")
  cat("  Selected items: ", summary_scalar(nrow(values$selected_items %||% data.frame())), "\n", sep = "")
  cat("  Exogenous scope: all fitted\n")
  cat("  Score total: all fitted items\n")
  cat("  Score cuts: ", cm2_cm3_score_cuts_label(values$score_cuts %||% integer()), "\n", sep = "")
  cat("  Two-way margins: ", summary_scalar(values$n_two_way_margins %||% 0L), "\n", sep = "")
  if (isTRUE(include_cm3)) {
    cat("  Three-way margins: ", summary_scalar(values$n_three_way_margins %||% 0L), "\n", sep = "")
  }

  cm2_row <- if (isTRUE(include_cm3)) values$cm2 else values$aggregate
  print_cm2_cm3_stat_block("CM2", cm2_row)
  if (isTRUE(include_cm3)) {
    print_cm2_cm3_stat_block("CM3", values$cm3)
  }
  print_cm2_cm3_stat_block("Item-trait", values$item_trait)
  cat("  Invariance terms: ", summary_scalar(nrow(values$invariance %||% data.frame())), "\n\n", sep = "")
  print_cm2_cm3_bootstrap_status(values$bootstrap)
  if (isTRUE(include_cm3)) {
    cat("Use summary(x) to show two-way margins, three-way margins, and decompositions.\n")
  } else {
    cat("Use summary(x) to show two-way margins and invariance by exogenous variable.\n")
  }
  invisible(x)
}

#' Internal print cm2 cm3 stat block helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param label Internal `label` value used by this helper.
#' @param row Internal `row` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_cm2_cm3_stat_block <- function(label, row) {
  row <- row %||% data.frame()
  if (!is.data.frame(row) || nrow(row) == 0L) {
    chi <- df <- p <- NA
  } else {
    chi <- row$chi_square[[1L]]
    df <- row$degrees_of_freedom[[1L]]
    p <- row$p_value[[1L]]
  }
  cat("\n")
  cat("  ", label, " Chisq: ", cm2_cm3_summary_chisq(chi), "\n", sep = "")
  cat("  ", label, " Df: ", summary_scalar(df), "\n", sep = "")
  cat("  ", label, " Pr(>Chisq): ", summary_p_value(p), "\n", sep = "")
  if (is.data.frame(row) && nrow(row) && "bootstrap_p_value" %in% names(row)) {
    cat(
      "  ", label, " Bootstrap Pr: ",
      summary_p_value(row$bootstrap_p_value[[1L]], empirical = TRUE),
      "\n",
      sep = ""
    )
  }
}

#' Internal print cm2 cm3 bootstrap status helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param bootstrap Whether to run bootstrap calibration.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_cm2_cm3_bootstrap_status <- function(bootstrap) {
  if (!is.list(bootstrap) || !isTRUE(bootstrap$enabled)) {
    cat("  Parametric bootstrap: not requested\n\n")
    return(invisible(NULL))
  }
  if (!isTRUE(bootstrap$possible)) {
    cat("  Parametric bootstrap: not possible\n")
    if (length(bootstrap$reasons)) {
      cat("  Reason: ", paste(bootstrap$reasons, collapse = "; "), "\n", sep = "")
    }
    cat("\n")
    return(invisible(NULL))
  }
  cat(
    "  Parametric bootstrap: ", bootstrap$nused, "/", bootstrap$nsim,
    " samples accepted (delta < ", summary_scalar(bootstrap$acceptance_delta), ")\n\n",
    sep = ""
  )
  invisible(NULL)
}

#' Internal cm2 cm3 summary chisq helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `cm2_cm3_summary_chisq()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_summary_chisq <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  format(signif(x, dig_tst), digits = dig_tst)
}

#' @export
summary.gRm_cm2 <- function(object, ...) {
  reject_summary_which(...)
  new_cm2_cm3_summary(object, include_cm3 = FALSE)
}

#' @export
summary.gRm_cm3 <- function(object, ...) {
  reject_summary_which(...)
  new_cm2_cm3_summary(object, include_cm3 = TRUE)
}

#' Internal new cm2 cm3 summary helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @param include_cm3 Internal `include_cm3` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_cm2_cm3_summary <- function(object, include_cm3) {
  values <- object$values %||% list()
  title <- if (isTRUE(include_cm3)) "gRm: CM3 fit diagnostic" else "gRm: CM2 fit diagnostic"
  tables <- list(
    aggregates = cm2_cm3_aggregate_table(values, include_cm3 = include_cm3),
    item_trait = cm2_cm3_item_trait_table(values),
    invariance = cm2_cm3_invariance_table(values),
    margins = cm2_cm3_margin_summary_table(values, include_cm3 = include_cm3)
  )
  bh <- rbind_fill(
    values$cm2_bh %||% data.frame(),
    if (isTRUE(include_cm3)) values$cm3_bh %||% data.frame() else data.frame()
  )
  if (is.data.frame(bh) && nrow(bh)) {
    tables$bh <- data.frame(
      Diagnostic = bh$diagnostic,
      FDR = bh$fdr,
      `Critical p` = bh$critical_p,
      `P-values` = as.integer(bh$n_p_values),
      `Source limit reached` = bh$source_limit_reached,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (isTRUE(values$bootstrap$enabled)) {
    tables$bootstrap <- cm2_cm3_bootstrap_summary_table(values$bootstrap)
  }
  new_gRm_summary(
    object,
    title = title,
    which = names(tables),
    tables = tables,
    header = c(
      paste0(
        "Selection: ", cm2_cm3_selection_metadata(values)$mode,
        "; focal items: ", paste(values$selected_items$item_name %||% character(), collapse = ", "),
        "; exogenous: all fitted; score total: all fitted items"
      ),
      paste0("Score cuts: ", cm2_cm3_score_cuts_label(values$score_cuts %||% integer()))
    ),
    table_names = names(tables),
    print_table_names = FALSE
  )
}

#' Internal CM2/CM3 selection metadata compatibility helper
#'
#' New results store the authoritative selection record in `values$selection`.
#' The fallback keeps print and summary methods usable for older serialized or
#' manually constructed result objects.
#' @param values CM2/CM3 values list.
#' @return A structured selection metadata list.
#' @keywords internal
#' @noRd
cm2_cm3_selection_metadata <- function(values) {
  if (is.list(values$selection) && length(values$selection$mode) == 1L) {
    return(values$selection)
  }
  list(
    schema_version = 0L,
    mode = "legacy_unspecified",
    resolved_items = values$selected_items %||% data.frame(),
    exogenous_scope = "all_fitted",
    score_group_included = TRUE,
    score_total_scope = "all_fitted_items"
  )
}

#' Internal cm2 cm3 selection label helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `cm2_cm3_selection_label()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_selection_label <- function(values) {
  selected <- values$selected_items %||% data.frame()
  items <- selected$item_name %||% character()
  exogenous <- values$exogenous_names %||% character()
  paste(
    if (length(items)) paste(items, collapse = ", ") else "<none>",
    if (length(exogenous)) paste(exogenous, collapse = ", ") else "<none>",
    sep = " | "
  )
}

#' Internal cm2 cm3 score cuts label helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param score_cuts Resolved total-score cut values.
#' @return The internal `cm2_cm3_score_cuts_label()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_score_cuts_label <- function(score_cuts) {
  if (!length(score_cuts)) {
    return("<none>")
  }
  paste(as.integer(score_cuts), collapse = ", ")
}

#' Internal cm2 cm3 aggregate table helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @param include_cm3 Internal `include_cm3` value used by this helper.
#' @return The internal `cm2_cm3_aggregate_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_aggregate_table <- function(values, include_cm3) {
  rows <- if (isTRUE(include_cm3)) {
    rbind_fill(values$cm2 %||% data.frame(), values$cm3 %||% data.frame())
  } else {
    values$aggregate %||% data.frame()
  }
  cm2_cm3_public_stat_table(rows, first_column = "Diagnostic")
}

#' Internal cm2 cm3 item trait table helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `cm2_cm3_item_trait_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_item_trait_table <- function(values) {
  cm2_cm3_public_stat_table(values$item_trait %||% data.frame(), first_column = "Diagnostic")
}

#' Internal cm2 cm3 invariance table helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `cm2_cm3_invariance_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_invariance_table <- function(values) {
  table <- values$invariance %||% data.frame()
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    Exogenous = table$background_name,
    Chisq = table$chi_square,
    Df = as.integer(table$degrees_of_freedom),
    `Pr(>Chisq)` = table$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if ("bootstrap_p_value" %in% names(table)) {
    out$`Bootstrap Pr` <- table$bootstrap_p_value
    out$`Bootstrap n` <- as.integer(table$bootstrap_nused)
  }
  out
}

#' Internal cm2 cm3 margin summary table helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @param include_cm3 Internal `include_cm3` value used by this helper.
#' @return The internal `cm2_cm3_margin_summary_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_summary_table <- function(values, include_cm3) {
  tests <- values$tests %||% data.frame()
  if (!isTRUE(include_cm3) && is.data.frame(tests) && nrow(tests)) {
    tests <- tests[tests$is_cm2, , drop = FALSE]
  }
  if (!is.data.frame(tests) || !nrow(tests)) {
    return(data.frame())
  }
  out <- data.frame(
    Margin = tests$margin,
    Type = tests$margin_type,
    Item = cm2_cm3_blank_missing(tests$variable1),
    `Item/Exogenous` = cm2_cm3_blank_missing(tests$variable2),
    `Item/Exogenous/Score` = cm2_cm3_blank_missing(tests$variable3),
    Chisq = tests$chi_square,
    Df = as.integer(tests$degrees_of_freedom),
    `Pr(>Chisq)` = tests$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if ("bootstrap_p_value" %in% names(tests)) {
    out$`Bootstrap Pr` <- tests$bootstrap_p_value
    out$`Bootstrap n` <- as.integer(tests$bootstrap_nused)
  }
  out
}

#' Internal cm2 cm3 blank missing helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `cm2_cm3_blank_missing()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_blank_missing <- function(x) {
  x <- as.character(x %||% character())
  x[is.na(x)] <- ""
  x
}

#' Internal cm2 cm3 public stat table helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param first_column Internal `first_column` value used by this helper.
#' @return The internal `cm2_cm3_public_stat_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_public_stat_table <- function(table, first_column) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    Diagnostic = table$diagnostic,
    Chisq = table$chi_square,
    Df = as.integer(table$degrees_of_freedom),
    `Pr(>Chisq)` = table$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if ("bootstrap_p_value" %in% names(table)) {
    out$`Bootstrap Pr` <- table$bootstrap_p_value
    out$`Bootstrap n` <- as.integer(table$bootstrap_nused)
  }
  names(out)[[1L]] <- first_column
  out
}

#' Internal cm2 cm3 bootstrap summary table helper
#'
#' Supports the api cm2 cm3 implementation while preserving its internal contract.
#' @param bootstrap Whether to run bootstrap calibration.
#' @return The internal `cm2_cm3_bootstrap_summary_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_bootstrap_summary_table <- function(bootstrap) {
  if (!is.list(bootstrap) || !isTRUE(bootstrap$enabled)) {
    return(data.frame())
  }
  data.frame(
    Status = bootstrap$source_status %||% NA_character_,
    Possible = bootstrap$possible %||% NA,
    Simulations = as.integer(bootstrap$nsim %||% 0L),
    Accepted = as.integer(bootstrap$nused %||% 0L),
    # RandSeed is a uint32 bit pattern and may exceed R's signed-integer range.
    Seed = as.numeric(bootstrap$seed %||% NA_real_),
    Reestimated = bootstrap$reestimate %||% NA,
    Workers = as.integer(bootstrap$bootstrap_jobs %||% 0L),
    `Acceptance delta` = bootstrap$acceptance_delta %||% 0.1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
