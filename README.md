# shapCopula

**Semiparametric Inference for Conditional Shapley Feature Importance**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20392903.svg)](https://doi.org/10.5281/zenodo.20392903)
[![arXiv](https://img.shields.io/badge/arXiv-2609.10313-b31b1b.svg)](https://doi.org/10.48550/arXiv.2609.10313)


R implementation accompanying:

> Gnasso, A. (2026). *Semiparametric Inference for Conditional Shapley Feature Importance*. arXiv:2609.10313. <https://doi.org/10.48550/arXiv.2609.10313>

---

## What this package does

`shapCopula` computes **Conditional SAGE** (Shapley Additive Global importancE) values for a black-box predictor $f$, together with asymptotically valid Wald confidence intervals and hypothesis tests.

The conditional value function $v(S;\, x) = E[f(X) \mid X_S = x_S]$ is estimated by a **Gaussian copula with empirical marginals** (default) or a **nonparametric vine copula** (`rvinecopulib`). A **one-step debiased estimator with K-fold cross-fitting** and a U-statistic debiasing of the squared loss yield root-n consistent, asymptotically normal estimates under mild regularity conditions (Theorem 1 in the paper).

A **marginal (interventional) SAGE** baseline is included for comparison with standard TreeSHAP-based attributions.

---

## Installation

```r
# Install from GitHub
remotes::install_github("agostinognasso/shapCopula")
```

**Dependencies:** `rvinecopulib >= 0.6.0`, R >= 4.1.0.

```r
install.packages("rvinecopulib")
```

---

## Quick start

```r
library(shapCopula)

# n x p data matrix + black-box predictor
f_predict <- function(Xm) predict(my_model, Xm)

# Conditional SAGE with 95% CIs (Gaussian copula, K=3 folds)
result <- sage_cv_gcop(f_predict, X, K = 3, M = 30, B = 20, seed = 42)
print(result)
#   feature  psi_hat      se   ci_lo   ci_hi      z  pvalue
#        X1   0.8312  0.0421  0.7487  0.9137  19.74 < 2e-16
#        X2   0.2145  0.0318  0.1522  0.2768   6.74   2e-11
#        X3   0.0031  0.0105 -0.0175  0.0237   0.30   0.768
#       ...
```

Columns returned:

| Column | Description |
|--------|-------------|
| `psi_hat` | SAGE estimate $\hat\Psi_j$ |
| `se` | Asymptotic standard error |
| `ci_lo`, `ci_hi` | $(1-\alpha)$ Wald confidence interval |
| `z` | Wald statistic for $H_0: \Psi_j = 0$ |
| `pvalue` | Two-sided p-value |

---

## Main functions

| Function | Sampler | Scalability |
|----------|---------|-------------|
| `sage_cv_gcop()` | Gaussian copula + empirical marginals | $p \lesssim 15$, **recommended** |
| `sage_cv()` | Nonparametric vine copula (TLL) | $p \lesssim 10$ |
| `sage_cv_marginal()` | Marginal bootstrap (interventional) | any $p$, baseline |

Supporting functions: `fit_gauss_copula()`, `fit_copula()`.

---

## Reproducing paper experiments

The simulation and application scripts live in the parent research repository,
under `reproduction/` and `Applicazione/`.

```bash
Rscript reproduction/run_simA_calibration.R              # calibration
Rscript reproduction/run_simB_power.R                    # power
Rscript reproduction/run_simC_correlation.R              # coverage vs correlation
Rscript reproduction/run_simD_misspecification.R         # copula misspecification
Rscript reproduction/run_simE_conditional_vs_marginal.R  # conditional vs marginal
Rscript Applicazione/run_applications_v2.R               # UCI Concrete, California Housing
```

Default tuning: $K = 3$, $M = 30$, $B = 20$. Sensitivity is studied in Appendix C.
One replication takes about 13 s at $n = 500$, $p = 5$; set `R_REPS` to shorten a run.

---

## Computational notes

- Complexity is $O(K \times M \times p \times n \times B)$.
- Gaussian copula fitting is $O(np^2)$ — much cheaper than vine copula fitting, which scales as $O(p^2)$ pair-copulas.
- The conditional sampler uses two independent sub-batches of size $B/2$ per step (U-statistic debiasing) to remove the Monte Carlo squared-noise bias.
- For $p > 15$, runtime becomes substantial. The package does not currently parallelize across features; use `parallel::mclapply` externally if needed.

---

## For developers

```r
# Regenerate NAMESPACE and man/ from roxygen2 comments
devtools::document()

# Run tests
devtools::test()

# Build vignette
devtools::build_vignettes()
```

---

## Citation

```bibtex
@unpublished{Gnasso2026condSAGE,
  author = {Agostino Gnasso},
  title  = {Inference for Conditional {Shapley} Values via Vine Copulas},
  year   = {2026},
  note   = {Working Paper}
}
```

For the software specifically:

```bibtex
@software{shapCopula,
  author  = {Agostino Gnasso},
  title   = {{shapCopula}: Semiparametric Inference for Conditional Shapley Feature Importance},
  year    = {2026},
  url     = {https://github.com/agostinognasso/shapCopula}
}
```
