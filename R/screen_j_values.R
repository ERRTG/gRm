#' Derive DIGRAM SCREEN J item-screening values
#'
#' Computes the source-shaped values printed by DIGRAM's `SCREEN J` command.
#' This is the item-screening branch from `DIGRAM1f.pas` command 3 and
#' `SKbias7.pas::Item_screening`; it is distinct from the graphical
#' `SCREEN` command, which is outside the current package scope.
#'
#' The validation runtime target uses the asymptotic, two-sided branch by default.
#' Exact mode follows the source `GENTAB1` Monte Carlo branch and computes
#' partial p-values directly from the project data.
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param project A parsed DIGRAM project.
#' @param exact Logical; when `TRUE`, replace the partial item-item and partial
#'   DIF p-values with the source-shaped Monte Carlo exact p-values used by
#'   `SKbias3.XYZ_bias_ANALYSE`.
#' @param nsim Number of random tables for exact p-values.
#' @param seed Random seed for the source-shaped exact branch.
#' @return A `gRm_screen_j_values` object.
#' @param repeated Whether to use repeated sequential simulation.
#' @param exact_state Internal `exact_state` value used by this helper.
#' @keywords internal
#' @noRd
screen_j_values <- function(project,
                            exact = FALSE,
                            repeated = FALSE,
                            nsim = 1000L,
                            seed = 9L,
                            exact_state = NULL) {
  exact <- isTRUE(exact)
  repeated <- isTRUE(repeated)
  nsim <- as.integer(nsim)
  seed <- as.integer(seed)
  if (is.null(exact_state)) {
    exact_state <- gRm_exact_state_from_flags(exact, repeated, nsim = nsim, seed = seed)
  }
  exact <- isTRUE(exact_state$exact)
  repeated <- isTRUE(exact_state$sequential)
  nsim <- as.integer(exact_state$nsim)
  seed <- as.integer(exact_state$seed)
  if (exact && (length(nsim) != 1L || is.na(nsim) || nsim < 1L)) {
    stop("`nsim` must be a positive integer for exact SCREEN J.", call. = FALSE)
  }
  if (exact && (length(seed) != 1L || is.na(seed))) {
    stop("`seed` must be a single integer for exact SCREEN J.", call. = FALSE)
  }
  seq_limit <- screen_j_repeated_seq_limit(project, repeated, nsim, exact_state = exact_state)
  seq_p0 <- if (isTRUE(repeated)) as.numeric(exact_state$seq_p0[[1L]]) else 0.05
  seq_boundary <- if (isTRUE(repeated)) as.numeric(exact_state$seq_b[[1L]]) else 1.058
  old_seed <- NULL
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (exact) {
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  items <- project$items
  backgrounds <- project$backgrounds
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  largest_possible_score <- sum(items$raw_max - 1L)

  item_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = n_items)
  for (item_index in seq_len(n_items)) {
    item_matrix[, item_index] <- as.integer(project$raw_data[, items$position[[item_index]]])
  }
  complete_items <- rep(TRUE, nrow(item_matrix))
  for (item_index in seq_len(n_items)) {
    complete_items <- complete_items &
      item_matrix[, item_index] >= 1L &
      item_matrix[, item_index] <= items$raw_max[[item_index]]
  }
  item_score <- rowSums(sweep(item_matrix, 2L, 1L, "-"))
  max_score <- if (any(complete_items)) max(item_score[complete_items]) else largest_possible_score
  good_items <- complete_items & item_score >= 1L &
    item_score <= max_score &
    item_score <= (largest_possible_score - 1L)
  item_screen <- screen_j_item_screening_values(
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    good_items = good_items,
    items = items,
    max_score = max_score,
    exact = exact,
    repeated = repeated,
    nsim = nsim,
    seed = seed,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )

  attr(project, "screen_j_exact") <- exact
  attr(project, "screen_j_nsim") <- nsim
  attr(project, "screen_j_seed") <- seed
  attr(project, "screen_j_repeated") <- repeated
  attr(project, "screen_j_seq_limit") <- seq_limit
  attr(project, "screen_j_seq_p0") <- seq_p0
  attr(project, "screen_j_seq_boundary") <- seq_boundary
  attr(project, "screen_j_exact_state") <- exact_state

  marginal_exo <- screen_j_exo_values(
    project = project,
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    items = items,
    backgrounds = backgrounds,
    max_score = max_score,
    partial = FALSE
  )
  partial_exo <- screen_j_exo_values(
    project = project,
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    items = items,
    backgrounds = backgrounds,
    max_score = max_score,
    partial = TRUE
  )

  bh_item_values <- c(as.vector(item_screen$partial_p[
    row(item_screen$partial_p) != col(item_screen$partial_p)
  ]))
  bh_values <- bh_item_values
  if (n_exo > 0L) {
    bh_values <- c(bh_values, as.vector(partial_exo$p))
  }
  weighted_gamma <- screen_j_weighted_partial_gamma(
    item_screen$partial_ppq,
    item_screen$partial_pmq
  )
  initial_item_bias <- partial_exo$p <= source_bh_critical(bh_values, 0.05)
  post_screen_dif <- screen_j_post_screen_dif(
    item_matrix = item_matrix,
    exo_values = partial_exo$exo_values,
    item_score = item_score,
    complete_items = complete_items,
    items = items,
    backgrounds = backgrounds,
    # Source trace: SKbias13.AnalysisOfSpuriousItemBias and
    # StepwiseItemBiasAnalysis dimension their score variable by
    # Largest_possible_score, not by the largest observed score.
    max_score = largest_possible_score,
    partial_exo = partial_exo,
    initial_item_bias = initial_item_bias,
    fdr_01 = source_bh_critical(bh_values, 0.01),
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
  local_dependence <- screen_j_stepwise_local_dependence(
    item_screen$partial_gamma,
    item_screen$partial_p,
    weighted_gamma,
    source_bh_critical(bh_values, 0.05)
  )
  score_effects <- screen_j_score_effect_values(
    project = project,
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    backgrounds = backgrounds,
    max_score = max_score,
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
  screen_j_assemble_values(
    items = items,
    backgrounds = backgrounds,
    item_screen = item_screen,
    marginal_exo = marginal_exo,
    partial_exo = partial_exo,
    weighted_gamma = weighted_gamma,
    local_dependence = local_dependence,
    initial_item_bias = initial_item_bias,
    post_screen_dif = post_screen_dif,
    score_effects = score_effects,
    bh_values = bh_values,
    exact = exact,
    nsim = nsim,
    seed = seed
  )
}

#' Compute the item-by-item SCREEN J matrices
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' Mathematical step: traverse row item then column item in source order,
#' retaining the asymmetric complete/good-record policy and exact-test branch.
#' @param item_matrix One-based source item values.
#' @param item_score Zero-based total item scores.
#' @param complete_items Complete-item record mask.
#' @param good_items Source interior-score record mask.
#' @param items Parsed item metadata.
#' @param max_score Largest complete-record observed score.
#' @param exact Whether exact conditional inference is active.
#' @param repeated Whether repeated sequential inference is active.
#' @param nsim Maximum number of simulated tables.
#' @param seed Source RNG seed.
#' @param seq_limit Sequential simulation limit.
#' @param seq_p0 Sequential null boundary.
#' @param seq_boundary Sequential stopping boundary.
#' @return Named marginal, partial, and item-restscore statistic arrays.
#' @keywords internal
#' @noRd
screen_j_item_screening_values <- function(item_matrix,
                                           item_score,
                                           complete_items,
                                           good_items,
                                           items,
                                           max_score,
                                           exact,
                                           repeated,
                                           nsim,
                                           seed,
                                           seq_limit,
                                           seq_p0,
                                           seq_boundary) {
  n_items <- nrow(items)
  marginal_gamma <- matrix(0, nrow = n_items, ncol = n_items)
  marginal_p <- matrix(2, nrow = n_items, ncol = n_items)
  partial_gamma <- matrix(0, nrow = n_items, ncol = n_items)
  partial_p <- matrix(2, nrow = n_items, ncol = n_items)
  partial_ppq <- matrix(0, nrow = n_items, ncol = n_items)
  partial_pmq <- matrix(0, nrow = n_items, ncol = n_items)
  rest_gamma <- rep(0, n_items)
  rest_p <- rep(2, n_items)

  dimnames(marginal_gamma) <- dimnames(marginal_p) <-
    dimnames(partial_gamma) <- dimnames(partial_p) <-
    dimnames(partial_ppq) <- dimnames(partial_pmq) <-
    list(items$label_code, items$label_code)
  names(rest_gamma) <- names(rest_p) <- items$label_code

  for (row_item in seq_len(n_items)) {
    # Source trace: SKbias7.inexpensive_rosenberg sets THE_ITEMS(.X.) := 1
    # before calling SCORE(THE_ITEMS), so the score column is Score \ X, not
    # Score minus the raw item code.
    rest_score <- item_score - (item_matrix[, row_item] - 1L)
    rest_table <- screen_j_pair_table(
      item_matrix[, row_item],
      rest_score + 1L,
      items$raw_max[[row_item]],
      max_score + 1L,
      complete_items
    )
    rest_stats <- screen_rc_gamma(rest_table)
    rest_gamma[[row_item]] <- rest_stats$gamma
    rest_p[[row_item]] <- rest_stats$p_value

    for (col_item in seq_len(n_items)) {
      if (row_item == col_item) {
        next
      }
      marginal_table <- screen_j_pair_table(
        item_matrix[, row_item],
        item_matrix[, col_item],
        items$raw_max[[row_item]],
        items$raw_max[[col_item]],
        complete_items
      )
      marginal_stats <- screen_rc_gamma(marginal_table)
      marginal_gamma[row_item, col_item] <- marginal_stats$gamma
      marginal_p[row_item, col_item] <- marginal_stats$p_value

      if (row_item > col_item) {
        x_values <- item_matrix[, col_item]
        y_values <- item_matrix[, row_item]
        x_dim <- items$raw_max[[col_item]]
        y_dim <- items$raw_max[[row_item]]
        partial_valid <- good_items
      } else {
        x_values <- item_matrix[, row_item]
        y_values <- item_matrix[, col_item]
        x_dim <- items$raw_max[[row_item]]
        y_dim <- items$raw_max[[col_item]]
        partial_valid <- complete_items
      }
      strata <- screen_j_strata_table(
        x_values,
        y_values,
        rest_score + 1L,
        x_dim,
        y_dim,
        max_score + 1L,
        partial_valid
      )
      partial_stats <- screen_j_partial_gamma(strata)
      partial_gamma[row_item, col_item] <- partial_stats$gamma
      partial_p[row_item, col_item] <- if (exact) {
        native_pair <- if (screen_j_item_pair_native_allowed(repeated, seq_p0, seq_boundary)) {
          screen_j_item_pair_conditional_exact_native(
            x_values,
            y_values,
            x_dim,
            y_dim,
            rest_score + 1L,
            max_score + 1L,
            partial_valid,
            nsim,
            seed,
            sequential = repeated,
            seq_limit = seq_limit
          )
        } else {
          NULL
        }
        if (!is.null(native_pair)) {
          native_pair$p_gamma_exact
        } else {
          partial_chi <- screen_j_partial_chi(strata)
          screen_j_exact_chi_gamma_prepared_r(
            screen_j_prepare_exact_slices(screen_j_strata_slices(strata)),
            partial_chi$stat,
            partial_stats$gamma,
            nsim,
            seed,
            sequential = repeated,
            seq_limit = seq_limit,
            seq_p0 = seq_p0,
            seq_boundary = seq_boundary
          )$p_gamma
        }
      } else {
        partial_stats$p_value
      }
      partial_ppq[row_item, col_item] <- partial_stats$ppq
      partial_pmq[row_item, col_item] <- partial_stats$pmq
    }
  }

  list(
    marginal_gamma = marginal_gamma,
    marginal_p = marginal_p,
    partial_gamma = partial_gamma,
    partial_p = partial_p,
    partial_ppq = partial_ppq,
    partial_pmq = partial_pmq,
    rest_gamma = rest_gamma,
    rest_p = rest_p
  )
}

#' Assemble a SCREEN J value object
#'
#' Keeps public object construction separate from the source-order numerical
#' traversal in `screen_j_item_screening_values()`.
#' @param items Parsed item metadata.
#' @param backgrounds Parsed background metadata.
#' @param item_screen Item-pair screening arrays.
#' @param marginal_exo Marginal item-by-background statistics.
#' @param partial_exo Partial item-by-background statistics.
#' @param weighted_gamma Weighted partial-gamma matrix.
#' @param local_dependence Selected local-dependence results.
#' @param initial_item_bias Initial DIF selection matrix.
#' @param post_screen_dif Post-screen DIF analysis.
#' @param score_effects Score-effect analysis.
#' @param bh_values Combined source BH p-value vector.
#' @param exact Whether exact inference was used.
#' @param nsim Requested simulation count.
#' @param seed Exact-test seed.
#' @return A `gRm_screen_j_values` object.
#' @keywords internal
#' @noRd
screen_j_assemble_values <- function(items,
                                     backgrounds,
                                     item_screen,
                                     marginal_exo,
                                     partial_exo,
                                     weighted_gamma,
                                     local_dependence,
                                     initial_item_bias,
                                     post_screen_dif,
                                     score_effects,
                                     bh_values,
                                     exact,
                                     nsim,
                                     seed) {
  structure(
    list(
      items = items,
      backgrounds = backgrounds,
      marginal = list(
        item_gamma = item_screen$marginal_gamma,
        item_p = item_screen$marginal_p,
        rest_gamma = item_screen$rest_gamma,
        rest_p = item_screen$rest_p,
        exo_stat = marginal_exo$stat,
        exo_p = marginal_exo$p,
        exo_kind = marginal_exo$kind,
        score_stat = marginal_exo$score_stat,
        score_p = marginal_exo$score_p
      ),
      partial = list(
        item_gamma = item_screen$partial_gamma,
        item_p = item_screen$partial_p,
        item_ppq = item_screen$partial_ppq,
        item_pmq = item_screen$partial_pmq,
        weighted_gamma = weighted_gamma,
        # Source trace: SKbias7.Item_screening prints
        # sum(abs(part_g[i,j] + part_g[j,i])) / (nitems * (nitems - 1)).
        average_abs_gamma = screen_j_average_abs_partial_gamma(item_screen$partial_gamma),
        exo_stat = partial_exo$stat,
        exo_p = partial_exo$p,
        exo_ppq = partial_exo$ppq,
        exo_pmq = partial_exo$pmq,
        exo_kind = partial_exo$kind
      ),
      model = list(
        local_dependence = local_dependence,
        initial_item_bias = initial_item_bias,
        item_bias = post_screen_dif$item_bias,
        item_bias_status = post_screen_dif$item_bias_status,
        spurious_dif = post_screen_dif$spurious_dif,
        multiple_dif = post_screen_dif$multiple_dif,
        score_effects = score_effects
      ),
      bh = list(
        fdr_05 = source_bh_critical(bh_values, 0.05),
        fdr_01 = source_bh_critical(bh_values, 0.01),
        fdr_001 = source_bh_critical(bh_values, 0.001),
        n_tests = length(bh_values)
      ),
      source_status = c(
        command = "DIGRAM1f_command_3_SCREEN_parameter_J",
        values = if (exact) {
          "SKbias7_Item_screening_inexpensive_rosenberg_inexpensive_itembias1_exact_GENTAB1"
        } else {
          "SKbias7_Item_screening_inexpensive_rosenberg_inexpensive_itembias1_asymptotic"
        },
        rendering = "SKbias7_Item_screening_WRITE_MESSAGE_cases_1_2_3"
      ),
      exact = exact,
      nsim = if (exact) nsim else 0L,
      seed = if (exact) seed else NA_integer_
    ),
    class = "gRm_screen_j_values"
  )
}
