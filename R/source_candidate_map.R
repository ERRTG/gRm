#' Internal source candidate map helper
#'
#' Supports the source candidate map implementation while preserving its internal contract.
#' @param n_candidates Internal `n_candidates` value used by this helper.
#' @param jobs Requested worker count.
#' @param fit_one Internal `fit_one` value used by this helper.
#' @return The internal `source_candidate_map()` computation result.
#' @keywords internal
#' @noRd
source_candidate_map <- function(n_candidates, jobs, fit_one) {
  n_candidates <- as.integer(n_candidates[[1L]])
  if (n_candidates <= 0L) {
    return(NULL)
  }
  jobs <- max(1L, min(as.integer(jobs[[1L]]), 128L, n_candidates))
  indices <- seq_len(n_candidates)
  rows <- if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(indices, fit_one, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    lapply(indices, fit_one)
  }
  do.call(rbind, rows)
}
