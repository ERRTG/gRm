lib_dir <- if (nzchar(R_ARCH)) {
  file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
} else {
  file.path(R_PACKAGE_DIR, "libs")
}
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

shared_objects <- Sys.glob(paste0("*", SHLIB_EXT))
if (length(shared_objects)) {
  copied <- file.copy(shared_objects, lib_dir, overwrite = TRUE)
  if (!all(copied)) {
    stop("failed to install shared library", call. = FALSE)
  }
  if (!WINDOWS) {
    Sys.chmod(file.path(lib_dir, basename(shared_objects)), mode = "755")
  }
}

if (file.exists("symbols.rds")) {
  file.copy("symbols.rds", lib_dir, overwrite = TRUE)
}

is_pkgload_compile <- grepl(
  "devtools_install_",
  normalizePath(R_PACKAGE_DIR, winslash = "/", mustWork = FALSE),
  fixed = TRUE
)
if (!is_pkgload_compile) {
  unlink(c(Sys.glob("*.o"), shared_objects, "symbols.rds"))
}
