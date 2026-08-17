#' Source-faithful DIGRAM exact/asymptotic command state
#'
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param inference One of asymptotic, exact, repeated, or sequential.
#' @param nsim Number of random tables for exact modes.
#' @param seed Random seed.
#' @param critlevel DIGRAM command critical level on the per-1000 scale.
#' @param risk DIGRAM repeated-Monte-Carlo risk on the per-1000 scale.
#' @return A `gRm_exact_command_state` list.
#' @keywords internal
#' @noRd
gRm_exact_command_state <- function(inference = c("asymptotic", "exact", "repeated", "sequential"),
                                       nsim = 1000L,
                                       seed = 9L,
                                       critlevel = 50L,
                                       risk = 1L) {
  critlevel_supplied <- !missing(critlevel)
  inference <- match.arg(inference)
  nsim <- exact_integer_scalar(nsim, "nsim", nonnegative = TRUE)
  seed <- exact_integer_scalar(seed, "seed")
  critlevel <- exact_integer_scalar(critlevel, "critlevel", nonnegative = TRUE)
  risk <- exact_integer_scalar(risk, "risk", nonnegative = TRUE)

  exact_nsim <- if (nsim == 0L) 1000L else nsim
  source_p0 <- source_seq_p0(critlevel)
  source_alpha <- source_seq_alpha(risk)
  state <- switch(
    inference,
    asymptotic = list(
      inference = "asymptotic",
      command = "ASYMPTOTIC",
      command_no = 17L,
      exact = FALSE,
      sequential = FALSE,
      nsim = 0L,
      seed = 9L,
      seq_p0 = NA_real_,
      seq_alpha = NA_real_,
      seq_b = NA_real_,
      seq_limit = 0L
    ),
    exact = list(
      inference = "exact",
      command = "EXA",
      command_no = 2L,
      exact = TRUE,
      sequential = FALSE,
      nsim = exact_nsim,
      seed = seed,
      seq_p0 = NA_real_,
      seq_alpha = NA_real_,
      seq_b = NA_real_,
      seq_limit = exact_nsim
    ),
    repeated = list(
      inference = "repeated",
      command = "REPEATED",
      command_no = 74L,
      exact = TRUE,
      sequential = TRUE,
      nsim = exact_nsim,
      seed = seed,
      seq_p0 = source_p0,
      seq_alpha = source_alpha,
      seq_b = source_seq_boundary(critlevel, risk, exact_nsim),
      # Source trace: PAS_scd/DIGRAM1f.pas command 74 (REP) sets the
      # repeated-Monte-Carlo boundary and keeps the same count cutoff as SEQ.
      seq_limit = 20L
    ),
    sequential = list(
      inference = "sequential",
      command = "SEQUENTIAL",
      command_no = 57L,
      exact = TRUE,
      sequential = TRUE,
      nsim = exact_nsim,
      seed = seed,
      seq_p0 = source_p0,
      seq_alpha = 0,
      seq_b = 1000,
      seq_limit = if (critlevel_supplied) critlevel else 20L
    )
  )

  class(state) <- c("gRm_exact_command_state", "list")
  state
}

#' Internal exact integer scalar helper
#'
#' Supports the exact command state implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param value Value to validate or transform.
#' @param name Internal name or label.
#' @param nonnegative Internal `nonnegative` value used by this helper.
#' @return The internal `exact_integer_scalar()` computation result.
#' @keywords internal
#' @noRd
exact_integer_scalar <- function(value, name, nonnegative = FALSE) {
  qualifier <- if (isTRUE(nonnegative)) {
    "single non-negative integer-like value"
  } else {
    "single integer-like value"
  }
  fail <- function() {
    stop("`", name, "` must be a ", qualifier, ".", call. = FALSE)
  }

  if (length(value) != 1L || !is.numeric(value)) {
    fail()
  }
  if (is.na(value) || !is.finite(value) || value != floor(value)) {
    fail()
  }
  if (isTRUE(nonnegative) && value < 0) {
    fail()
  }
  if (value < -.Machine$integer.max || value > .Machine$integer.max) {
    fail()
  }

  as.integer(value)
}

#' Internal gRm exact state from flags helper
#'
#' Supports the exact command state implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param repeated Whether to use repeated sequential simulation.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @return The internal `gRm_exact_state_from_flags()` computation result.
#' @keywords internal
#' @noRd
gRm_exact_state_from_flags <- function(exact, repeated, nsim = 1000L, seed = 9L) {
  inference <- if (!isTRUE(exact)) {
    "asymptotic"
  } else if (isTRUE(repeated)) {
    "repeated"
  } else {
    "exact"
  }
  gRm_exact_command_state(inference, nsim = nsim, seed = seed)
}

#' Internal gRm exact command state public helper
#'
#' Supports the exact command state implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param inference Internal `inference` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param critlevel Internal `critlevel` value used by this helper.
#' @param risk Internal `risk` value used by this helper.
#' @return The internal `gRm_exact_command_state_public()` computation result.
#' @keywords internal
#' @noRd
gRm_exact_command_state_public <- function(inference,
                                              nsim = 1000L,
                                              seed = 9L,
                                              critlevel = NULL,
                                              risk = NULL) {
  inference <- match.arg(inference, c("asymptotic", "exact", "repeated", "sequential"))
  if (identical(inference, "repeated")) {
    return(gRm_exact_command_state(
      "repeated",
      nsim = nsim,
      seed = seed,
      critlevel = if (is.null(critlevel)) 50L else critlevel,
      risk = if (is.null(risk)) 1L else risk
    ))
  }
  if (identical(inference, "sequential")) {
    if (is.null(critlevel)) {
      return(gRm_exact_command_state("sequential", nsim = nsim, seed = seed))
    }
    return(gRm_exact_command_state("sequential", nsim = nsim, seed = seed, critlevel = critlevel))
  }
  gRm_exact_command_state(inference, nsim = nsim, seed = seed)
}

#' Internal source seq p0 helper
#'
#' Supports the exact command state implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param critlevel Internal `critlevel` value used by this helper.
#' @return The internal `source_seq_p0()` computation result.
#' @keywords internal
#' @noRd
source_seq_p0 <- function(critlevel = 50L) {
  p0 <- exact_integer_scalar(critlevel, "critlevel", nonnegative = TRUE) / 1000
  if (p0 < 0.10) {
    0.05
  } else if (p0 < 0.25) {
    0.10
  } else {
    0.25
  }
}

#' Internal source seq alpha helper
#'
#' Supports the exact command state implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param risk Internal `risk` value used by this helper.
#' @return The internal `source_seq_alpha()` computation result.
#' @keywords internal
#' @noRd
source_seq_alpha <- function(risk = 1L) {
  alpha <- exact_integer_scalar(risk, "risk", nonnegative = TRUE) / 1000
  if (alpha < 0.005) {
    0.001
  } else if (alpha < 0.010) {
    0.005
  } else if (alpha < 0.020) {
    0.010
  } else if (alpha < 0.050) {
    0.020
  } else {
    0.05
  }
}

#' Source boundary used by DIGRAM's repeated Monte Carlo command
#'
#' Mirrors `SKrandom.pas::SEQ_INIT` for command 74. The source first snaps
#' `P0`, `ALPHA`, and `NSIM` to coarse table coordinates, then selects `SEQ_B`.
#'
#' Source trace: `source/PAS_skunits/SKrandom.pas::GENTAB1`.
#' @param critlevel Internal `critlevel` value used by this helper.
#' @param risk Internal `risk` value used by this helper.
#' @param nsim Requested simulation count.
#' @return The internal `source_seq_boundary()` computation result.
#' @keywords internal
#' @noRd
source_seq_boundary <- function(critlevel = 50L, risk = 1L, nsim = 1000L) {
  p0 <- source_seq_p0(critlevel)
  alpha <- source_seq_alpha(risk)
  nsim <- exact_integer_scalar(nsim, "nsim", nonnegative = TRUE)
  n <- (nsim %/% 100L) * 100L
  if (nsim - n > 0L) n <- n + 100L
  if (n > 1000L) n <- 1000L

  key <- paste(p0, alpha, sep = ":")
  table <- list(
    "0.05:0.001" = c(default = 1.058),
    "0.05:0.005" = c("100" = 0.849, "200" = 0.870, "300" = 0.872, "400" = 0.874, "500" = 0.875, "600" = 0.882, "700" = 0.883, "800" = 0.883, "900" = 0.883, "1000" = 0.883, default = 0.883),
    "0.05:0.01" = c("100" = 0.777, "200" = 0.810, "300" = 0.815, "400" = 0.823, "500" = 0.827, "600" = 0.828, "700" = 0.828, "800" = 0.828, "900" = 0.828, "1000" = 0.829, default = 0.883),
    "0.05:0.02" = c("100" = 0.702, "200" = 0.725, "300" = 0.739, "400" = 0.739, "500" = 0.743, "600" = 0.744, "700" = 0.749, "800" = 0.749, "900" = 0.750, "1000" = 0.752, default = 0.752),
    "0.05:0.05" = c("100" = 0.573, "200" = 0.603, "300" = 0.621, "400" = 0.633, "500" = 0.639, "600" = 0.642, "700" = 0.645, "800" = 0.645, "900" = 0.645, "1000" = 0.646, default = 0.646),
    "0.1:0.001" = c("100" = 1.287, "200" = 1.313, "300" = 1.317, "400" = 1.317, "500" = 1.317, "600" = 1.317, "700" = 1.321, "800" = 1.347, "900" = 1.347, "1000" = 1.347, default = 1.347),
    "0.1:0.005" = c("100" = 1.095, "200" = 1.116, "300" = 1.137, "400" = 1.147, "500" = 1.147, "600" = 1.162, "700" = 1.163, "800" = 1.163, "900" = 1.163, "1000" = 1.163, default = 1.163),
    "0.1:0.01" = c("100" = 1.000, "200" = 1.025, "300" = 1.042, "400" = 1.053, "500" = 1.061, "600" = 1.067, "700" = 1.067, "800" = 1.068, "900" = 1.068, "1000" = 1.068, default = 1.068),
    "0.1:0.02" = c("100" = 0.904, "200" = 0.941, "300" = 0.956, "400" = 0.971, "500" = 0.980, "600" = 0.983, "700" = 0.983, "800" = 0.990, "900" = 0.991, "1000" = 0.992, default = 0.992),
    "0.1:0.05" = c("100" = 0.772, "200" = 0.811, "300" = 0.828, "400" = 0.837, "500" = 0.850, "600" = 0.852, "700" = 0.855, "800" = 0.858, "900" = 0.861, "1000" = 0.863, default = 0.863),
    "0.25:0.001" = c("100" = 1.641, "200" = 1.686, "300" = 1.714, "400" = 1.721, "500" = 1.724, "600" = 1.724, "700" = 1.724, "800" = 1.724, "900" = 1.724, "1000" = 1.724, default = 1.724),
    "0.25:0.005" = c("100" = 1.471, "200" = 1.526, "300" = 1.539, "400" = 1.541, "500" = 1.552, "600" = 1.563, "700" = 1.563, "800" = 1.565, "900" = 1.565, "1000" = 1.566, default = 1.566),
    "0.25:0.01" = c("100" = 1.383, "200" = 1.427, "300" = 1.443, "400" = 1.461, "500" = 1.466, "600" = 1.470, "700" = 1.473, "800" = 1.477, "900" = 1.485, "1000" = 1.488, default = 1.488),
    "0.25:0.02" = c("100" = 1.257, "200" = 1.304, "300" = 1.337, "400" = 1.350, "500" = 1.354, "600" = 1.362, "700" = 1.371, "800" = 1.371, "900" = 1.376, "1000" = 1.379, default = 1.379),
    "0.25:0.05" = c("100" = 1.066, "200" = 1.133, "300" = 1.157, "400" = 1.174, "500" = 1.186, "600" = 1.192, "700" = 1.200, "800" = 1.205, "900" = 1.207, "1000" = 1.210, default = 1.210)
  )
  row <- table[[key]]
  if (is.null(row)) return(2)
  value <- row[as.character(n)]
  if (is.na(value)) value <- row["default"]
  as.numeric(value)
}
