# Source-shaped Goodman-Kruskal/RC gamma table statistics.
#
# Source trace: source/GLLRM.txt::PREPARE_REAL_GAMMA_STATISTICS uses the
# all-cell convention. For each table cell it builds AIJ and DIJ by scanning
# every other table cell, then accumulates P = sum(tab * AIJ),
# Q = sum(tab * DIJ), PPQ = P + Q, and PMQ = P - Q.
source_rc_gamma_stats <- function(tab, include_cells = FALSE) {
  tab <- as.matrix(tab)
  storage.mode(tab) <- "double"

  if (length(dim(tab)) != 2L || nrow(tab) < 1L || ncol(tab) < 1L) {
    stop("`tab` must be a non-empty two-way table.", call. = FALSE)
  }

  aij <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
  dij <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
  p <- 0
  q <- 0

  for (row in seq_len(nrow(tab))) {
    for (col in seq_len(ncol(tab))) {
      for (other_row in seq_len(nrow(tab))) {
        for (other_col in seq_len(ncol(tab))) {
          concordant <- (row > other_row && col > other_col) ||
            (row < other_row && col < other_col)
          discordant <- (row < other_row && col > other_col) ||
            (row > other_row && col < other_col)
          if (concordant) {
            aij[row, col] <- aij[row, col] + tab[other_row, other_col]
          } else if (discordant) {
            dij[row, col] <- dij[row, col] + tab[other_row, other_col]
          }
        }
      }
      p <- p + tab[row, col] * aij[row, col]
      q <- q + tab[row, col] * dij[row, col]
    }
  }

  ppq <- p + q
  pmq <- p - q
  out <- list(
    p = p,
    q = q,
    ppq = ppq,
    pmq = pmq,
    gamma = if (ppq > 0) pmq / ppq else 0
  )

  if (include_cells) {
    out$aij <- aij
    out$dij <- dij
  }

  out
}

source_rc_gamma_counts <- function(tab) {
  stats <- source_rc_gamma_stats(tab, include_cells = FALSE)
  list(gamma = stats$gamma, ppq = stats$ppq, pmq = stats$pmq)
}

# Source trace: source/PAS_skunits/skfit2.pas::Standardize_tab4 performs
# 30 fixed row/column scaling passes. Pascal stores table cells as
# V1V2_TAB(.I,J.) and scales the first source table index to Cmarg before the
# second source table index to Rmarg. R callers pass those source indices as
# matrix rows and columns, so this helper scales rows before columns and keeps
# the fixed pass count instead of iterating to convergence.
source_standardize_table_margins <- function(tab, row_margins, col_margins, n_iter = 30L) {
  out <- as.matrix(tab)
  storage.mode(out) <- "double"
  row_margins <- as.numeric(row_margins)
  col_margins <- as.numeric(col_margins)

  if (length(row_margins) != nrow(out) || length(col_margins) != ncol(out)) {
    stop("Margin lengths must match the table dimensions.", call. = FALSE)
  }

  for (step in seq_len(n_iter)) {
    row_sums <- rowSums(out)
    for (row in seq_along(row_sums)) {
      if (row_sums[[row]] > 0) {
        out[row, ] <- out[row, ] * row_margins[[row]] / row_sums[[row]]
      }
    }

    col_sums <- colSums(out)
    for (col in seq_along(col_sums)) {
      if (col_sums[[col]] > 0) {
        out[, col] <- out[, col] * col_margins[[col]] / col_sums[[col]]
      }
    }
  }

  out
}

gamma_cell_tables <- function(tab) {
  stats <- source_rc_gamma_stats(tab, include_cells = TRUE)
  list(aij = stats$aij, dij = stats$dij, p = stats$p, q = stats$q)
}
