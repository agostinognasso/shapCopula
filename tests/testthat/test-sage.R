set.seed(1)
n <- 200; p <- 4
X <- matrix(stats::rnorm(n * p), n, p)
colnames(X) <- paste0("X", seq_len(p))
# Linear model: only X1 and X2 matter
f_lin <- function(Xm) 3 * Xm[, 1] - 2 * Xm[, 2]

# ---- fit_gauss_copula -------------------------------------------------------

test_that("fit_gauss_copula returns expected structure", {
  gc <- fit_gauss_copula(X)
  expect_named(gc, c("ecdfs", "sorted_cols", "Sigma", "X_train"))
  expect_equal(dim(gc$Sigma), c(p, p))
  expect_true(all(diag(gc$Sigma) == 1))
  expect_true(isSymmetric(gc$Sigma))
  expect_length(gc$ecdfs, p)
  expect_length(gc$sorted_cols, p)
})

test_that("fit_gauss_copula Sigma is positive definite", {
  gc <- fit_gauss_copula(X)
  expect_true(all(eigen(gc$Sigma, only.values = TRUE)$values > 0))
})

# ---- sage_cv_gcop -----------------------------------------------------------

result_gcop <- sage_cv_gcop(f_lin, X, K = 2L, M = 15L, B = 8L,
                             alpha = 0.05, seed = 42)

test_that("sage_cv_gcop returns correct data.frame structure", {
  expect_s3_class(result_gcop, "data.frame")
  expect_equal(nrow(result_gcop), p)
  expect_named(result_gcop, c("feature", "psi_hat", "se", "ci_lo", "ci_hi",
                               "z", "pvalue"))
})

test_that("sage_cv_gcop CIs are well-formed", {
  expect_true(all(result_gcop$ci_lo < result_gcop$ci_hi))
  expect_true(all(result_gcop$se > 0))
  expect_true(all(result_gcop$pvalue >= 0 & result_gcop$pvalue <= 1))
})

test_that("sage_cv_gcop concentrates importance on the signal features", {
  # The *population* parameter satisfies Psi_j >= 0, since conditioning on more
  # features cannot increase the residual variance. The *estimator* carries no
  # such constraint: for X3 and X4 the truth is exactly zero, so the estimates
  # straddle zero and asserting non-negativity there would only be asserting
  # that noise happens to fall on one side. What must hold is that the signal
  # features are significantly positive and dominate the null ones.
  psi <- stats::setNames(result_gcop$psi_hat, result_gcop$feature)
  signal <- result_gcop$feature %in% c("X1", "X2")
  expect_true(all(psi[c("X1", "X2")] > 0))
  expect_true(all(result_gcop$ci_lo[signal] > 0))
  expect_lt(max(abs(psi[c("X3", "X4")])), min(psi[c("X1", "X2")]) / 2)
})

test_that("sage_cv_gcop detects X1 as more important than X3/X4", {
  # X1 has coefficient 3, X2 has -2; X3, X4 are noise
  ord <- order(result_gcop$psi_hat, decreasing = TRUE)
  # X1 (index 1) and X2 (index 2) should rank above X3 and X4
  top2 <- sort(result_gcop$feature[ord[1:2]])
  expect_true(all(c("X1", "X2") %in% top2))
})

# ---- sage_cv_marginal -------------------------------------------------------

result_marg <- sage_cv_marginal(f_lin, X, K = 2L, M = 15L, B = 8L,
                                alpha = 0.05, seed = 42)

test_that("sage_cv_marginal returns correct data.frame structure", {
  expect_s3_class(result_marg, "data.frame")
  expect_equal(nrow(result_marg), p)
  expect_named(result_marg, c("feature", "psi_hat", "se", "ci_lo", "ci_hi",
                               "z", "pvalue"))
})

test_that("sage_cv_marginal detects active features", {
  ord  <- order(result_marg$psi_hat, decreasing = TRUE)
  top2 <- sort(result_marg$feature[ord[1:2]])
  expect_true(all(c("X1", "X2") %in% top2))
})

# ---- Independence: conditonal ≈ marginal for iid X -------------------------

test_that("gcop and marginal SAGE agree on the active set for iid data", {
  # Under independence the conditional and the marginal value function coincide,
  # so both estimators must single out the same active features. Correlating the
  # full rank vectors is not a meaningful check here: X3 and X4 are both exactly
  # null, so their relative order is pure noise and the Spearman coefficient
  # moves with it.
  top2 <- function(res) sort(res$feature[order(res$psi_hat, decreasing = TRUE)][1:2])
  expect_equal(top2(result_gcop), c("X1", "X2"))
  expect_equal(top2(result_marg), top2(result_gcop))
})

# ---- fit_copula / sage_cv (vine sampler) ------------------------------------

set.seed(7)
Xv <- matrix(stats::rnorm(150 * 3), 150, 3)
colnames(Xv) <- paste0("X", seq_len(3))
f_v <- function(Xm) 2 * Xm[, 1] - Xm[, 2]

test_that("fit_copula returns a vinecop object alongside empirical marginals", {
  vc <- fit_copula(Xv, family_set = "parametric")
  expect_named(vc, c("vc", "ecdfs", "sorted_cols", "X_train", "n_train"))
  expect_s3_class(vc$vc, "vinecop")
  expect_length(vc$ecdfs, 3)
  expect_length(vc$sorted_cols, 3)
  expect_identical(dim(vc$X_train), c(150L, 3L))
  expect_identical(vc$n_train, 150L)
})

test_that("fit_copula honours a fixed D-vine order", {
  vc <- fit_copula(Xv, family_set = "parametric", order = c(3L, 1L, 2L))
  expect_identical(as.integer(rvinecopulib::get_structure(vc$vc)$order),
                   c(3L, 1L, 2L))
})

test_that("sage_cv with the vine sampler recovers the active features", {
  res <- sage_cv(f_v, Xv, K = 2L, M = 8L, B = 6L,
                 family_set = "parametric", seed = 3)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("feature", "psi_hat", "se", "ci_lo", "ci_hi",
                      "z", "pvalue"))
  expect_equal(nrow(res), 3)
  expect_true(all(res$se > 0))
  expect_true(all(res$ci_lo < res$ci_hi))
  # f depends on X1 (coefficient 2) and X2 (-1); X3 is pure noise.
  expect_equal(res$feature[which.max(res$psi_hat)], "X1")
  expect_gt(res$ci_lo[res$feature == "X1"], 0)
})


# ---- S3 methods -------------------------------------------------------------

test_that("the estimators return a sage_estimate that is still a data.frame", {
  expect_s3_class(result_gcop, "sage_estimate")
  expect_s3_class(result_gcop, "data.frame")
  expect_identical(attr(result_gcop, "sampler"), "Gaussian copula")
  expect_identical(attr(result_gcop, "K"), 2L)
  # subsetting and coercion must keep working for downstream code
  expect_type(result_gcop$psi_hat, "double")
  expect_equal(nrow(as.data.frame(result_gcop)), p)
})

test_that("print and plot methods run without error", {
  expect_output(print(result_gcop), "Conditional SAGE")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(result_gcop))
})

test_that("print and plot accept the arguments they forward to `...`", {
  # Regression: both methods used to fix row.names / xlim / ylim / ylab and
  # forward `...` at the same time, so supplying any of them raised
  # "formal argument matched by multiple actual arguments".
  expect_output(print(result_gcop, row.names = TRUE), "Conditional SAGE")
  expect_output(print(result_gcop, row.names = FALSE), "Conditional SAGE")
  expect_output(print(result_gcop, digits = 5), "Conditional SAGE")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(result_gcop, xlim = c(-1, 1)))
  expect_silent(plot(result_gcop, ylim = c(0, p + 1)))
  expect_silent(plot(result_gcop, ylab = "feature"))
  expect_silent(plot(result_gcop, xlim = c(-1, 1), ylab = "feature"))
})
