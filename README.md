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

# All items (default), asymptotic
cm2(fit0)
cm3(fit0)

# Selected items, asymptotic
cm2(fit0, items = c("I1", "I3"))
cm3(fit0, items = c("I1", "I3"))

# All items (default), parametric bootstrap
cm2(fit0, bootstrap = TRUE, nsim = 100, seed = 47)
cm3(fit0, bootstrap = TRUE, nsim = 100, seed = 47)

# Selected items, parametric bootstrap
cm2(fit0, items = c("I1", "I3"), bootstrap = TRUE, nsim = 100, seed = 47)
cm3(fit0, items = c("I1", "I3"), bootstrap = TRUE, nsim = 100, seed = 47)
```

The installed package computes numeric results from data and model objects
through native R functions.

## Public API

The exported functions are:

- `gRm()` and `read_digram_project()` for data setup;
- `gllrm()` and `fit()` for model specification and estimation;
- `model_graph()` for graph extraction from model and fit objects;
- `screen()` and `score_effects()` for screening and exogenous score-effect
  diagnostics;
- `item_fit()`, `local_dependence()`, `dif()`, `global_homogeneity()`, and
  `ari()` for post-fit numeric results;
- source-backed `cm2()` and `cm3()` fit diagnostics, with opt-in parametric
  bootstrap calibration.

For CM2/CM3, `items = NULL` selects all fitted items, matching DIGRAM's blank
item prompt. Exact item names or one-based indices select at least two focal
items and are normalized to fitted source order. This selection changes only
the reported diagnostic margins: all fitted exogenous variables remain
automatic, score groups use the total across all fitted items, fitted LD/DIF
terms retain unselected partners, and bootstrap generation/refitting continues
to use the complete fitted model. The supplied fit is not mutated.

Fitted item parameters and thresholds are reported through `summary(fit)` and
`summary(fit, which = "parameters" / "thresholds")`. Some diagnostics, such as
`item_fit()` and `score_effects()`, return their public tables directly. Use
`summary()` for analysis, model, fit, and diagnostic result objects, and use
`which =` only for documented multi-view outputs, such as fitted parameters,
item-fit tables, and global-homogeneity sections.

## Source Faithfulness

The R implementation is a parallel implementation of selected non-GUI DIGRAM
computations in native R. Exported R functions and production computational
helpers compute results directly from the original source material.

The current package version covers ordinal item and exogenous variables. Items
and exogenous variables are interpreted as ordinal source variables.

Legacy DIGRAM import support covers the simple ordinal category coding subset:
category codes in `DIGRAM.imv` must be contiguous one-based codes matching the
values in `DIGRAM.csv`.

Where a historical DIGRAM value lacks source-backed provenance, R results use
`NA` plus source-status metadata. This includes the global-homogeneity residual
cells whose source formula is not available.

## Development Notes

Run the package tests with:

```sh
R CMD INSTALL .
Rscript -e "testthat::test_local(reporter = 'summary')"
```

Production R code implements algorithms directly from the original source
material.
