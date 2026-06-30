list_to_one_row <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (!is.list(x) || !length(x)) {
    return(data.frame())
  }
  data.frame(as.list(x), check.names = FALSE, stringsAsFactors = FALSE)
}

result_matrix_table <- function(x, key_name, value_name) {
  if (!is.matrix(x) || !length(x)) {
    return(data.frame())
  }
  idx <- which(!is.na(x), arr.ind = TRUE)
  out <- data.frame(
    item = rownames(x)[idx[, 1L]],
    key = colnames(x)[idx[, 2L]],
    value = as.numeric(x[idx]),
    stringsAsFactors = FALSE
  )
  stats::setNames(out, c("item", key_name, value_name))
}

result_named_numeric_table <- function(x, name_col, value_col) {
  if (!length(x)) {
    return(data.frame())
  }
  out <- data.frame(
    name = names(x),
    value = as.numeric(x),
    stringsAsFactors = FALSE
  )
  stats::setNames(out, c(name_col, value_col))
}

rbind_fill_typed_na <- function(x) {
  if (is.integer(x)) {
    return(NA_integer_)
  }
  if (is.numeric(x)) {
    return(NA_real_)
  }
  if (is.logical(x)) {
    return(NA)
  }
  if (is.character(x)) {
    return(NA_character_)
  }
  x[NA_integer_]
}

rbind_fill <- function(...) {
  parts <- list(...)
  parts <- parts[vapply(parts, nrow, integer(1L)) > 0L]
  if (!length(parts)) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  missing_values <- stats::setNames(lapply(cols, function(col) {
    for (part in parts) {
      if (col %in% names(part)) {
        return(rbind_fill_typed_na(part[[col]]))
      }
    }
    NA
  }), cols)
  parts <- lapply(parts, function(x) {
    missing <- setdiff(cols, names(x))
    for (col in missing) {
      x[[col]] <- missing_values[[col]]
    }
    x[cols]
  })
  do.call(rbind, parts)
}

normalize_summary_table <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) {
    return(data.frame())
  }
  rownames(x) <- NULL
  x
}
