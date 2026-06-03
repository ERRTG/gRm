# gRm

`gRm` is a native R implementation of the GLLRM-oriented parts of DIGRAM.
The public API is a compact statistical modeling workflow:

```r
analysis <- gRm(data, items = c("I1", "I2", "I3"), exogenous = "group")
model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:group)
fit0 <- fit(model)

summary(fit0)
summary(item_parameters(fit0), which = "coefficients")
summary(item_fit(fit0), which = "tests")
summary(local_dependence(fit0), which = "tests")
summary(dif(fit0), which = "tests")
summary(global_homogeneity(fit0), which = "tests")
summary(score_effects(analysis), which = "tests")
```

The installed package computes numeric results from data and model objects. It
does not generate DIGRAM fixed-width output files, parse historical DIGRAM
runtime output, or compare package output to oracle files. Repository-level
validation against historical DIGRAM output lives under
`validation/digram_oracle/`.

## Public API

The exported functions are:

- `gRm()` and `read_digram_project()` for data setup;
- `gllrm()` and `fit()` for model specification and estimation;
- `screen()` and `score_effects()` for screening and exogenous score-effect
  diagnostics;
- `item_parameters()`, `item_fit()`, `local_dependence()`, `dif()`, and
  `global_homogeneity()` for post-fit numeric results.

`summary()` is the public presentation surface. Use `which =` to select a
named summary section, such as `"fit"`, `"coefficients"`, `"tests"`, or
`"bh"`, depending on the object.

## Source Faithfulness

The R implementation is a parallel implementation of non-GUI DIGRAM, not a
wrapper around the Pascal harness. Exported R functions and production
computational helpers must compute results in R and must not invoke Pascal
binaries, Pascal harness scripts, shell commands, generated Pascal TSV outputs,
or cached Pascal results.

When a historical DIGRAM value is not source-backed, R exposes `NA` plus
unmodeled metadata rather than fitting constants from oracle output. In
particular, global-homogeneity residual cells without available source
provenance remain unmodeled.

## Development Notes

Run the package tests with:

```sh
R CMD INSTALL gRm
Rscript -e "library(gRm); testthat::test_dir('gRm/tests/testthat', reporter='summary')"
```

Production R code must keep implementing algorithms directly from the original
source material. Oracle report parsing and semantic comparison belong in
`validation/digram_oracle/`, not in the installed package namespace.
