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

test_that("sage_cv_gcop SAGE values sum is non-negative for squared loss", {
  # SAGE values (squared-loss based) are non-negative by construction
  expect_true(all(result_gcop$psi_hat >= -1e-3))
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

test_that("gcop and marginal SAGE agree qualitatively on iid data", {
  # Feature rankings should be consistent (Spearman rank correlation > 0.5)
  r <- stats::cor(result_gcop$psi_hat, result_marg$psi_hat,
                  method = "spearman")
  expect_gt(r, 0.5)
})
