# shapCopula 0.2.1

## Fixes

* `print()` and `plot()` on a `sage_estimate` fixed `row.names`, and `xlim`,
  `ylim` and `ylab` respectively, while also forwarding `...` to the underlying
  method. Supplying any of those arguments, all of them documented as passed
  through `...`, raised "formal argument matched by multiple actual arguments".
  They are now treated as defaults that `...` overrides, and the test suite
  covers each of them.

# shapCopula 0.2.0

## Conditional sampling

* `sage_cv()` now draws conditional samples through the Rosenblatt transform of
  a D-vine. One D-vine is fitted per Shapley permutation, with an order chosen
  so that every coalition visited along that permutation is a suffix of the
  order and is therefore exactly conditionable. This replaces the earlier
  importance-sampling scheme, which reweighted a pool of joint draws with a
  Gaussian kernel, and it changes the values `sage_cv()` returns.
* `fit_copula()` gained an `order` argument that fixes the D-vine structure,
  which is what makes the above possible.

## New methods

* `print()` and `plot()` methods for the `sage_estimate` class returned by the
  estimators. `print()` shows the Wald summary under a header giving the sampler
  and the tuning parameters; `plot()` draws a forest plot, one interval per
  feature, with features whose interval excludes zero drawn solid. The object
  remains a plain data frame underneath.

## Fixes

* `testthat.R` sat in `tests/testthat/` rather than in `tests/`, so `R CMD check`
  never ran the test suite. It is in the right place now, and the tests have
  been extended to cover the exact vine sampler, the S3 methods and the marginal
  baseline.

## Other

* `inst/CITATION` added, pointing at the accompanying preprint.
* README rewritten around a worked case study on the UCI Concrete data.
* Package title aligned with the accompanying paper.

# shapCopula 0.1.0

First release.

* `sage_cv_gcop()` estimates conditional SAGE values with a Gaussian copula and
  empirical marginals, returning Wald confidence intervals and two-sided tests
  from a one-step cross-fitted estimator.
* `sage_cv()` does the same with a vine copula, drawing conditional samples by
  importance reweighting of a pool of joint draws.
* `sage_cv_marginal()` provides the marginal (interventional) baseline, which
  coincides with the conditional estimator under independence and diverges from
  it under dependence.
* `fit_gauss_copula()` and `fit_copula()` expose the underlying copula fits.
