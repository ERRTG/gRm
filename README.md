# gRm

`gRm` is a native R implementation of the GLLRM-oriented parts of DIGRAM.
The public API is a compact statistical modeling workflow:

```r
analysis <- gRm(data, items = c("I1", "I2", "I3"), exogenous = "site")
model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site)
fit0 <- fit(model)

fit0
summary(fit0)
summary(fit0, which = "parameters")
summary(fit0, which = "thresholds")

item_fit(fit0)
item_fit(fit0, which = "items")

ld <- local_dependence(fit0)
ld
summary(ld)

dif_tests <- dif(fit0)
dif_tests
summary(dif_tests)

global <- global_homogeneity(fit0)
global
summary(global)

score_effects(analysis)

m2(fit0) # experimental
m3(fit0) # experimental
```

The installed package computes numeric results from data and model objects. It
does not generate DIGRAM fixed-width output files, parse historical DIGRAM
runtime output, or compare package output to oracle files. Repository-level
validation against historical DIGRAM output lives under
`validation/digram_oracle/`.

## Public API

The exported functions are:

- `gRm()` and `read_digram_project()` for data setup;
- `gllrm()`, `model_graph()`, `plot()`, `fit()`, `logLik()`, and `anova()`
  for model specification, graph inspection, estimation, and likelihood
  comparison;
- `screen()` and `score_effects()` for screening and exogenous score-effect
  diagnostics;
- `item_fit()`, `local_dependence()`, `dif()`, `global_homogeneity()`, and
  `ari()` for post-fit numeric results;
- experimental `m2()` and `m3()` fit diagnostics.

Fitted item parameters and thresholds are reported through `summary(fit)` and
`summary(fit, which = "parameters" / "thresholds")`; there is no separate
`item_parameters()` accessor. Some diagnostics, such as `item_fit()` and
`score_effects()`, return their public tables directly. Use `summary()` for
analysis, model, fit, and diagnostic result objects, and use `which =` only for
documented multi-view outputs, such as fitted parameters, item-fit tables, and
global-homogeneity sections.

## Source Faithfulness

The R implementation is a parallel implementation of non-GUI DIGRAM, not a
wrapper around the Pascal harness. Exported R functions and production
computational helpers must compute results in R and must not invoke Pascal
binaries, Pascal harness scripts, shell commands, generated Pascal TSV outputs,
or cached Pascal results.

The current package version is ordinal-only. Items and exogenous variables are
interpreted as ordinal source variables. Nominal and mixed DIGRAM variable-type
behavior is not implemented, so multi-category nominal exogenous variables are
outside the currently source-faithful package scope.

Legacy DIGRAM import support is also restricted to the simple ordinal category
coding subset: category codes in `DIGRAM.imv` must be contiguous one-based codes
matching the values in `DIGRAM.csv`. Historical projects with zero-based,
non-contiguous, or separately recoded `.imv` category mappings are not
source-faithfully implemented in the current package.

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
