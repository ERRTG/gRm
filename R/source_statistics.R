#' Source normal tail probability
#'
#' Port of `SkStat.TAILNORM`, used by the source chi-square tail approximation.
#'
#' Source trace: `source/PAS_scd/SkStat.pas::TAILNORM`.
#' @param value Normal deviate.
#' @param upper Logical; `TRUE` returns the upper tail.
#' @return Tail probability.
#' @keywords internal
#' @noRd
source_tail_norm <- function(value, upper = TRUE) {
  # SkStat.TAILNORM uses these literal binary64 constants and changes the tail
  # flag, rather than the deviate, for a negative input.
  normal_c <- 0.3989422804014327
  big_x <- 170
  if (value == 0) {
    return(0.5)
  }
  if (value < 0) {
    upper <- !upper
  }
  value <- abs(value)
  x2 <- value * value
  if ((x2 / 2) < big_x) {
    y <- normal_c * exp(-0.5 * x2)
  } else {
    y <- 0
  }
  n <- y / value
  if (!upper && n == 0) {
    return(1)
  }
  if (upper && n == 0) {
    return(0)
  }
  if ((upper && value > 2.32) || (!upper && value > 3.5)) {
    # TAILNORM's large-deviate branch evaluates its two-term continued
    # fraction until both successive floating-point ratios stop changing.
    q1 <- value
    p2 <- y * value
    i <- 1
    p1 <- y
    q2 <- x2 + i
    if (upper) {
      s <- p1 / q1
      m <- s
      t <- p2 / q2
    } else {
      s <- i - p1 / q1
      m <- s
      t <- i - p2 / q2
    }
    while (m != t && s != t) {
      i <- i + 1
      s <- value * p2 + i * p1
      p1 <- p2
      p2 <- s
      s <- value * q2 + i * q1
      q1 <- q2
      q2 <- s
      s <- m
      m <- t
      if (upper) {
        t <- p2 / q2
      } else {
        t <- 1 - p2 / q2
      }
    }
    return(t)
  }

  s <- y * value
  term <- y * value
  i <- 1
  t <- 0
  while (s != t) {
    # The central branch is the source odd-power series; retain its exact
    # equality stop and addition order rather than calling stats::pnorm().
    i <- i + 2
    t <- s
    term <- term * x2 / i
    s <- s + term
  }
  if (upper) {
    0.5 - s
  } else {
    0.5 + s
  }
}

#' Source chi-square upper-tail probability
#'
#' Port of `SkStat.PFCHI`.
#'
#' Source trace: `source/PAS_scd/SkStat.pas::PFCHI`.
#' @param df Degrees of freedom.
#' @param x Chi-square statistic.
#' @return Upper-tail probability.
#' @keywords internal
#' @noRd
source_pfchi <- function(df, x) {
  source_co <- 0.3989422804014327
  big_x <- 170
  df_div_2 <- df %/% 2L
  if (df == 0L || x <= 0) {
    return(1)
  }
  if (df > 100L) {
    # SkStat.PFCHI uses its Wilson--Hilferty normal transform above 100 df and
    # passes that deviate through the source TAILNORM approximation.
    transformed <- sqrt(4.5 * df) * (exp(log(x / df) / 3) + 1 / (4.5 * df) - 1)
    return(source_tail_norm(transformed, TRUE))
  }
  if ((df %% 2L) == 0L) {
    # Even df: start at exp(-x/2), then add the source finite gamma-series
    # terms in increasing integer order.
    if (x < big_x) {
      p <- exp(-0.5 * x)
    } else {
      p <- 0
    }
    c_value <- p
    if (df_div_2 > 1L) {
      for (i in seq_len(df_div_2 - 1L)) {
        c_value <- c_value * 0.5 * x / i
        p <- p + c_value
      }
    }
    return(p)
  }

  p <- 2 * source_tail_norm(sqrt(x), TRUE)
  # Odd df: start with twice the source normal tail and add the half-integer
  # recurrence. This is deliberately not stats::pchisq().
  if (x < big_x) {
    c_value <- exp(-0.5 * x) * 2 * source_co / sqrt(x)
  } else {
    c_value <- 0
  }
  if (df_div_2 > 0L) {
    for (i in seq_len(df_div_2)) {
      c_value <- c_value * x / (2 * i - 1)
      p <- p + c_value
    }
  }
  p
}

#' Source Benjamini-Hochberg critical p-value
#'
#' Port of `SourceBenjaminiHochbergCritical`, matching the current Pascal
#' harness behavior.
#'
#' Source trace: `source/PAS_skunits/SkStat.pas::BenjaminiHochberg0`.
#' @param p_values Numeric p-values.
#' @param alpha False discovery rate target.
#' @return Critical p-value.
#' @keywords internal
#' @noRd
source_bh_critical <- function(p_values, alpha = 0.05) {
  n <- length(p_values)
  if (n <= 0L) {
    return(0)
  }
  sorted <- sort(p_values)
  result <- alpha / n
  # BenjaminiHochberg0 scans from the largest rank downward and assigns the
  # rank boundary before testing p_(i) <= i*alpha/n. If none passes, the last
  # assigned value is alpha/n; preserve that source boundary convention.
  for (index in rev(seq_len(n))) {
    candidate <- (index / n) * alpha
    result <- candidate
    if (sorted[[index]] <= candidate) {
      break
    }
  }
  result
}
