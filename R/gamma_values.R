#' Derive DIGRAM marginal gamma report values
#'
#' Computes the Goodman-Kruskal marginal association-gamma matrix printed by
#' DIGRAM's `GAMMA` command. The source path is
#' `DGRexe.execute_marginal_gamma`, which builds one two-variable marginal
#' hypothesis for every ordinal pair and stores `results[1,5]`. That value is
#' produced by `SkStat.RCGAMMA`; the same quadrant-count formula is mirrored in
#' `SourceRaschCore.SourceRCGammaStats`.
#'
#' Production R computes from the encoded `project$raw_data` matrix created by
#' either public constructor; Pascal and supplied DIGRAM reports are validation
#' oracles only.
#'
#' Source trace: `source/digram_source_20260817/skunits/SkStat.pas::RCGAMMA`.
#' @param project A parsed DIGRAM project from [read_digram_project()].
#' @return A `gRm_gamma_values` object.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' values <- gamma_values(project)
#' values$gamma[1:3, 1:3]
#' }
#' @keywords internal
#' @noRd
gamma_values <- function(project) {
  variables <- project$variables
  ordinal <- variables$raw_max > 1L & variables$vtype == 3L
  raw <- project$raw_data
  n_vars <- nrow(variables)
  gamma <- matrix(999, nrow = n_vars, ncol = n_vars)
  ppq <- matrix(0, nrow = n_vars, ncol = n_vars)
  pmq <- matrix(0, nrow = n_vars, ncol = n_vars)
  pair_n <- matrix(0L, nrow = n_vars, ncol = n_vars)
  dimnames(gamma) <- list(variables$label_code, variables$label_code)
  dimnames(ppq) <- dimnames(gamma)
  dimnames(pmq) <- dimnames(gamma)
  dimnames(pair_n) <- dimnames(gamma)
  # DIGRAM renders marginal gamma self-pairs as "." rather than as a numeric
  # association; expose the R-facing numeric value as missing/not applicable.
  diag(gamma) <- NA_real_

  for (row_index in seq_len(n_vars - 1L)) {
    for (col_index in seq.int(row_index + 1L, n_vars)) {
      if (!ordinal[[row_index]] || !ordinal[[col_index]]) {
        next
      }
      tab <- gRm_pairwise_gamma_table(
        raw[, variables$position[[row_index]]],
        raw[, variables$position[[col_index]]],
        variables$raw_max[[row_index]],
        variables$raw_max[[col_index]]
      )
      stats <- gRm_goodman_kruskal_gamma(tab)
      gamma[row_index, col_index] <- gamma[col_index, row_index] <- stats$gamma
      ppq[row_index, col_index] <- ppq[col_index, row_index] <- stats$ppq
      pmq[row_index, col_index] <- pmq[col_index, row_index] <- stats$pmq
      pair_n[row_index, col_index] <- pair_n[col_index, row_index] <- sum(tab)
    }
  }

  structure(
    list(
      variables = data.frame(
        label_code = variables$label_code,
        name = variables$name,
        raw_max = variables$raw_max,
        vtype = variables$vtype,
        ordinal = ordinal,
        stringsAsFactors = FALSE
      ),
      gamma = gamma,
      ppq = ppq,
      pmq = pmq,
      n = pair_n
    ),
    class = "gRm_gamma_values"
  )
}

#' Build a pairwise complete ordinal table for marginal gamma
#'
#' Source trace: `source/digram_source_20260817/skunits/SkStat.pas::RCGAMMA`.
#' @param row_values Integer source-coded values for the row variable.
#' @param col_values Integer source-coded values for the column variable.
#' @param row_dim Number of ordinal row categories.
#' @param col_dim Number of ordinal column categories.
#' @return Integer matrix of pairwise complete counts.
#' @keywords internal
#' @noRd
gRm_pairwise_gamma_table <- function(row_values, col_values, row_dim, col_dim) {
  valid <- row_values >= 1L & row_values <= row_dim & col_values >= 1L & col_values <= col_dim
  tab <- matrix(0L, nrow = row_dim, ncol = col_dim)
  for (index in which(valid)) {
    tab[row_values[[index]], col_values[[index]]] <- tab[row_values[[index]], col_values[[index]]] + 1L
  }
  tab
}

#' Goodman-Kruskal gamma statistics using DIGRAM's RC gamma formula
#'
#' Source trace: `source/digram_source_20260817/skunits/SkStat.pas::RCGAMMA`.
#' @param tab A two-way ordinal count table.
#' @return List with `gamma`, `ppq`, and `pmq`.
#' @keywords internal
#' @noRd
gRm_goodman_kruskal_gamma <- function(tab) {
  source_rc_gamma_counts(tab)
}
