test_that("source RC gamma stats match the Pascal all-cell convention", {
  pascal_reference <- function(tab) {
    tab <- as.matrix(tab)
    aij <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
    dij <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
    p <- 0
    q <- 0

    for (i in seq_len(nrow(tab))) {
      for (j in seq_len(ncol(tab))) {
        for (k in seq_len(nrow(tab))) {
          for (l in seq_len(ncol(tab))) {
            if ((i > k && j > l) || (i < k && j < l)) {
              aij[i, j] <- aij[i, j] + tab[k, l]
            } else if ((i < k && j > l) || (i > k && j < l)) {
              dij[i, j] <- dij[i, j] + tab[k, l]
            }
          }
        }
        p <- p + tab[i, j] * aij[i, j]
        q <- q + tab[i, j] * dij[i, j]
      }
    }

    ppq <- p + q
    pmq <- p - q
    list(
      aij = aij,
      dij = dij,
      p = p,
      q = q,
      ppq = ppq,
      pmq = pmq,
      gamma = if (ppq > 0) pmq / ppq else 0
    )
  }

  tables <- list(
    matrix(c(1, 2, 3, 4), nrow = 2L),
    matrix(c(0, 3, 1, 2, 4, 0, 5, 1, 2), nrow = 3L),
    matrix(c(5, 0, 0, 0, 0, 1), nrow = 2L),
    matrix(c(2, 0, 1), nrow = 1L),
    matrix(c(2, 0, 1), ncol = 1L),
    matrix(0, nrow = 3L, ncol = 2L)
  )

  for (tab in tables) {
    expected <- pascal_reference(tab)
    observed <- source_rc_gamma_stats(tab, include_cells = TRUE)
    expect_equal(observed$p, expected$p, tolerance = 0)
    expect_equal(observed$q, expected$q, tolerance = 0)
    expect_equal(observed$ppq, expected$ppq, tolerance = 0)
    expect_equal(observed$pmq, expected$pmq, tolerance = 0)
    expect_equal(observed$gamma, expected$gamma, tolerance = 0)
    expect_equal(observed$aij, expected$aij, tolerance = 0)
    expect_equal(observed$dij, expected$dij, tolerance = 0)
  }
})

test_that("source RC gamma compact counts match full stats", {
  tab <- matrix(c(0, 2, 5, 1, 3, 0, 4, 1, 2), nrow = 3L)
  stats <- source_rc_gamma_stats(tab, include_cells = TRUE)
  counts <- source_rc_gamma_counts(tab)

  expect_equal(stats$p, 26, tolerance = 0)
  expect_equal(stats$q, 134, tolerance = 0)
  expect_equal(stats$ppq, 160, tolerance = 0)
  expect_equal(stats$pmq, -108, tolerance = 0)
  expect_equal(stats$gamma, -0.675, tolerance = 0)
  expect_equal(counts$gamma, stats$gamma, tolerance = 0)
  expect_equal(counts$ppq, stats$ppq, tolerance = 0)
  expect_equal(counts$pmq, stats$pmq, tolerance = 0)
})

test_that("source table standardization follows Pascal Standardize_tab4 pass order", {
  pascal_standardize_tab4_reference <- function(tab, row_margins, col_margins, n_iter = 30L) {
    out <- tab
    for (step in seq_len(n_iter)) {
      for (row in seq_len(nrow(out))) {
        total <- sum(out[row, ])
        if (total > 0) {
          out[row, ] <- out[row, ] * row_margins[[row]] / total
        }
      }
      for (col in seq_len(ncol(out))) {
        total <- sum(out[, col])
        if (total > 0) {
          out[, col] <- out[, col] * col_margins[[col]] / total
        }
      }
    }
    out
  }

  tab <- matrix(
    c(
      4, 0, 3,
      2, 5, 7,
      1, 6, 8,
      0, 2, 9
    ),
    nrow = 3L
  )
  row_margins <- c(12, 0, 18)
  col_margins <- c(4, 10, 6, 10)

  observed <- source_standardize_table_margins(tab, row_margins, col_margins)
  expected <- pascal_standardize_tab4_reference(tab, row_margins, col_margins)

  expect_equal(dim(observed), c(3L, 4L))
  expect_equal(observed, expected, tolerance = 1e-14)
  expect_equal(observed[2L, ], rep(0, 4L), tolerance = 0)
})

test_that("source gamma wrappers preserve caller-specific contracts", {
  tab <- matrix(c(0, 2, 5, 1, 3, 0, 4, 1, 2), nrow = 3L)
  expected <- source_rc_gamma_stats(tab, include_cells = TRUE)

  local_stats <- local_independence_source_gamma_counts(tab)
  expect_equal(local_stats$gamma, expected$gamma, tolerance = 0)
  expect_equal(local_stats$ppq, expected$ppq, tolerance = 0)
  expect_equal(local_stats$pmq, expected$pmq, tolerance = 0)

  dif_stats <- dif_tests_source_gamma_counts(tab)
  expect_equal(dif_stats$gamma, expected$gamma, tolerance = 0)
  expect_equal(dif_stats$ppq, expected$ppq, tolerance = 0)
  expect_equal(dif_stats$pmq, expected$pmq, tolerance = 0)

  screen_counts <- screen_rc_gamma_counts(tab)
  expect_equal(screen_counts$gamma, expected$gamma, tolerance = 0)
  expect_equal(screen_counts$ppq, expected$ppq, tolerance = 0)
  expect_equal(screen_counts$pmq, expected$pmq, tolerance = 0)

  screen_stats <- screen_rc_gamma(tab)
  expect_equal(screen_stats$gamma, expected$gamma, tolerance = 0)
  expect_equal(screen_stats$ppq, expected$ppq, tolerance = 0)
  expect_equal(screen_stats$pmq, expected$pmq, tolerance = 0)

  marginal_stats <- gRm_goodman_kruskal_gamma(tab)
  expect_equal(marginal_stats$gamma, expected$gamma, tolerance = 0)
  expect_equal(marginal_stats$ppq, expected$ppq, tolerance = 0)
  expect_equal(marginal_stats$pmq, expected$pmq, tolerance = 0)

  expect_equal(goodman_kruskal_gamma(tab), expected$gamma, tolerance = 0)
  expect_equal(source_gamma_from_table(tab), expected$gamma, tolerance = 0)

  cells <- gamma_cell_tables(tab)
  expect_equal(cells$aij, expected$aij, tolerance = 0)
  expect_equal(cells$dij, expected$dij, tolerance = 0)
  expect_equal(cells$p, expected$p, tolerance = 0)
  expect_equal(cells$q, expected$q, tolerance = 0)

  fitted <- fitted_gamma_stats(tab)
  expect_equal(fitted$gamma, expected$gamma, tolerance = 0)
  expect_equal(fitted$variance, 0.0470455078125, tolerance = 1e-15)
})

test_that("marginal gamma diagonal is not an applicable self-pair value", {
  project <- list(
    variables = data.frame(
      label_code = c("a", "b"),
      name = c("item_a", "item_b"),
      raw_max = c(2L, 2L),
      vtype = c(3L, 3L),
      position = c(1L, 2L),
      stringsAsFactors = FALSE
    ),
    raw_data = matrix(c(1L, 1L, 1L, 2L, 2L, 2L), ncol = 2L, byrow = TRUE)
  )

  values <- gamma_values(project)

  expect_true(all(is.na(diag(values$gamma))))
  expect_equal(values$gamma[["a", "b"]], 1, tolerance = 0)
})
