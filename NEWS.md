# shapCopula 0.1.0

First release.

## Features

* `sage_cv_gcop()` estimates conditional SAGE values with a Gaussian copula and
  empirical marginals, returning Wald confidence intervals and two-sided tests
  from a one-step cross-fitted estimator.
* `sage_cv()` does the same with a vine copula. Conditional draws are exact:
  one D-vine is fitted per Shapley permutation, with an order chosen so that
  every coalition of that permutation is conditionable through the Rosenblatt
  transform.
* `sage_cv_marginal()` provides the marginal (interventional) baseline, which
  coincides with the conditional estimator under independence and diverges from
  it under dependence.
* `fit_gauss_copula()` and `fit_copula()` expose the underlying copula fits.
  `fit_copula()` accepts an `order` argument that fixes the D-vine structure.
