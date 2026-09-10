# Core SAGE estimation machinery.
#
# All three exported functions (sage_cv_gcop, sage_cv, sage_cv_marginal) share
# the same cross-fitting skeleton (.sage_crossfit) and the same debiased
# psi-local computation (.debiased_psi). They differ only in the conditional
# sampler passed to .debiased_psi.


# --- Internal: debiased psi-local matrix ------------------------------------

# Compute the N x p matrix of per-observation local SAGE contributions
# using a U-statistic debiasing of the squared loss (Covert et al. 2020).
#
# `sampler_factory` is called once per Shapley permutation and returns the
# sampler to use for that permutation, with signature
# sampler_fn(S, x_S_matrix, B_per_row) -> matrix. Most samplers ignore the
# permutation; the vine sampler uses it to fit a D-vine whose order makes every
# coalition of that permutation exactly conditionable (see
# vine_cond_sampler_exact).
.debiased_psi <- function(f_predict, X_eval, f_X_eval, sampler_factory, M, Bh) {
  N <- nrow(X_eval)
  p <- ncol(X_eval)
  psi <- matrix(0, N, p)

  delta_unb <- function(S, sampler_fn) {
    if (length(S) == p) return(rep(0, N))
    if (length(S) == 0) {
      # v_empty = E[f(X)] is constant; N*Bh unconditional samples suffice.
      sa <- sampler_fn(integer(0), matrix(0, 1, 0), N * Bh)
      sb <- sampler_fn(integer(0), matrix(0, 1, 0), N * Bh)
      va <- mean(f_predict(sa))
      vb <- mean(f_predict(sb))
      return((f_X_eval - va) * (f_X_eval - vb))
    }
    sa <- sampler_fn(S, X_eval[, S, drop = FALSE], Bh)
    sb <- sampler_fn(S, X_eval[, S, drop = FALSE], Bh)
    pa <- f_predict(sa); pb <- f_predict(sb)
    va <- rowMeans(matrix(pa, N, Bh, byrow = TRUE))
    vb <- rowMeans(matrix(pb, N, Bh, byrow = TRUE))
    (f_X_eval - va) * (f_X_eval - vb)
  }

  for (m in seq_len(M)) {
    perm       <- sample(seq_len(p))
    sampler_fn <- sampler_factory(perm)
    S      <- integer(0)
    d_prev <- delta_unb(S, sampler_fn)
    for (k in seq_len(p)) {
      j     <- perm[k]
      S_new <- c(S, j)
      d_new <- delta_unb(S_new, sampler_fn)
      psi[, j] <- psi[, j] + (d_prev - d_new) / M
      S      <- S_new
      d_prev <- d_new
    }
  }
  psi
}


# --- Internal: cross-fitting skeleton ----------------------------------------

# Runs K-fold cross-fitting: for each fold, calls fold_fn(train_idx, test_idx, X)
# and assembles the global psi_local matrix, then computes Wald inference.
.sage_crossfit <- function(f_predict, X, K, alpha, fold_fn, M = NA, sampler = NULL) {
  X   <- as.matrix(X)
  n   <- nrow(X); p <- ncol(X)
  nms <- colnames(X) %||% paste0("X", seq_len(p))

  folds     <- sample(rep(seq_len(K), length.out = n))
  psi_local <- matrix(NA_real_, n, p)
  colnames(psi_local) <- nms

  for (k in seq_len(K)) {
    test_idx  <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)
    psi_local[test_idx, ] <- fold_fn(train_idx, test_idx, X)
  }

  psi_hat <- colMeans(psi_local)
  se      <- apply(psi_local, 2, stats::sd) / sqrt(n)
  z       <- psi_hat / se
  ci_lo   <- psi_hat - stats::qnorm(1 - alpha / 2) * se
  ci_hi   <- psi_hat + stats::qnorm(1 - alpha / 2) * se
  pvalue  <- 2 * (1 - stats::pnorm(abs(z)))

  out <- data.frame(
    feature = nms,
    psi_hat = psi_hat,
    se      = se,
    ci_lo   = ci_lo,
    ci_hi   = ci_hi,
    z       = z,
    pvalue  = pvalue,
    row.names = NULL
  )
  # `sage_estimate` only adds print/plot methods; the object stays a data.frame
  # so that existing code that subsets or writes it out keeps working.
  class(out) <- c("sage_estimate", "data.frame")
  attr(out, "K") <- K
  attr(out, "M") <- M
  attr(out, "alpha") <- alpha
  out
}


# --- Exported: Gaussian-copula SAGE ------------------------------------------

#' Conditional SAGE via Gaussian Copula with Cross-Fitting
#'
#' Computes global feature importance (SAGE) using the conditional value
#' function \eqn{v(S; x) = E[f(X) \mid X_S = x_S]} estimated via a Gaussian
#' copula with empirical marginals. Inference relies on a one-step debiased
#' estimator with K-fold cross-fitting.
#'
#' @param f_predict Function mapping an n-by-p numeric matrix to a numeric
#'   vector of predictions of length n.
#' @param X Numeric matrix or data frame, shape n-by-p.
#' @param K Integer. Number of cross-fitting folds (default 3).
#' @param M Integer. Number of random permutations per fold for the
#'   permutation-sampling approximation of the Shapley weights (default 30).
#' @param B Integer. Conditional samples per evaluation point per permutation
#'   step (default 20). Internally split into two sub-batches of size `B/2`
#'   for U-statistic debiasing of the squared loss.
#' @param alpha Numeric. Significance level for Wald confidence intervals
#'   (default 0.05).
#' @param seed Optional integer for `set.seed()` before fold assignment.
#'
#' @return A `data.frame` with one row per feature and columns:
#' \describe{
#'   \item{`feature`}{Feature name (from `colnames(X)` or `"X1"`, ...).}
#'   \item{`psi_hat`}{One-step SAGE estimate \eqn{\hat\Psi_j}.}
#'   \item{`se`}{Asymptotic standard error \eqn{\hat\sigma_j / \sqrt{n}}.}
#'   \item{`ci_lo`, `ci_hi`}{Wald \eqn{(1-\alpha)} confidence interval.}
#'   \item{`z`}{Wald z-statistic for \eqn{H_0: \Psi_j = 0}.}
#'   \item{`pvalue`}{Two-sided p-value.}
#' }
#'
#' @details
#' The Gaussian copula is fitted on each training fold via [shapCopula::fit_gauss_copula()].
#' Conditional distributions \eqn{X_{-S} \mid X_S = x_S} are sampled
#' analytically on the latent normal scale, then back-transformed through
#' empirical marginal quantiles.
#'
#' **Computational complexity** scales as \eqn{O(K \times M \times p \times n \times B)},
#' feasible for \eqn{p \lesssim 15} on standard hardware. The Gaussian copula
#' has no pair-copula fitting overhead, making it substantially faster than
#' the full vine estimator [shapCopula::sage_cv()].
#'
#' **Debiasing.** Each squared-loss term
#' \eqn{\delta_S(x) = (f(x) - v_S(x))^2} is estimated by a U-statistic that
#' uses two independent sub-batches of size `B/2`, removing the
#' Monte-Carlo-squared noise bias from the plug-in squared estimator.
#'
#' @references
#' Gnasso, A. (2026). *Semiparametric Inference for Conditional Shapley Feature Importance*.
#' arXiv preprint, doi:10.48550/arXiv.2609.10313
#'
#' Covert, I., Lundberg, S., & Lee, S.-I. (2020). Understanding global feature
#' contributions with additive importance measures. *NeurIPS*, 33.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 400; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("X", seq_len(p))
#' # f depends mainly on X1, moderately on X2
#' f_lin <- function(Xm) 2 * Xm[, 1] - Xm[, 2]
#' result <- sage_cv_gcop(f_lin, X, K = 3, M = 20, B = 10, seed = 1)
#' print(result)
#' }
#'
#' @seealso [shapCopula::sage_cv()] for vine-copula conditional sampling;
#'   [shapCopula::sage_cv_marginal()] for the interventional baseline.
#'
#' @export
sage_cv_gcop <- function(f_predict, X, K = 3L, M = 30L, B = 20L,
                         alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Bh <- max(2L, as.integer(B %/% 2))

  fold_fn <- function(train_idx, test_idx, X) {
    gcop    <- fit_gauss_copula(X[train_idx, , drop = FALSE])
    X_test  <- X[test_idx, , drop = FALSE]
    f_X_test <- f_predict(X_test)
    sampler <- function(S, x_S_mat, B_per_row)
      gcop_cond_sampler(gcop, S, x_S_mat, B_per_row)
    .debiased_psi(f_predict, X_test, f_X_test, function(perm) sampler, M, Bh)
  }

  res <- .sage_crossfit(f_predict, X, K, alpha, fold_fn, M = M)
  attr(res, "sampler") <- "Gaussian copula"
  res
}


# --- Exported: Vine-copula SAGE ----------------------------------------------

#' Conditional SAGE via Nonparametric Vine Copula with Cross-Fitting
#'
#' Computes global feature importance (SAGE) using the conditional value
#' function estimated via a vine copula fitted with `rvinecopulib`. Conditional
#' draws are exact: one D-vine is fitted per Shapley permutation, whose order is
#' chosen so that every coalition visited along that permutation can be
#' conditioned on through the Rosenblatt transform. Inference follows the same
#' one-step cross-fitted Wald procedure as [shapCopula::sage_cv_gcop()].
#'
#' @inheritParams sage_cv_gcop
#' @param family_set Character. Copula family set for `rvinecopulib::vinecop`.
#'   Default `"nonparametric"` (TLL). Use `"parametric"` for faster fitting.
#' @param cores Integer. Parallel cores for pair-copula fitting.
#'
#' @return A `data.frame` with the same columns as [shapCopula::sage_cv_gcop()].
#'
#' @details
#' Vine copula fitting scales as \eqn{O(p^2)} pair-copulas. With nonparametric
#' families the estimator is tractable for \eqn{p \lesssim 10}; consider
#' `family_set = "parametric"` or [shapCopula::sage_cv_gcop()] for larger problems.
#'
#' **Conditional sampling.** The coalitions visited along a Shapley permutation
#' are its prefixes. Fitting a D-vine whose variable order is the reversed
#' permutation therefore makes each of them a *suffix* of the order, and for a
#' D-vine the Rosenblatt transform conditions sequentially from the end of the
#' order. Fixing the corresponding components of the transform and redrawing the
#' rest as independent uniforms yields exact draws from
#' \eqn{X_{-S} \mid X_S = x_S}, with no kernel smoothing and no importance
#' weights. The cost is `M` vine fits per fold instead of one; the benefit is
#' that the sampler is unbiased, whereas kernel reweighting shrinks conditional
#' means toward the centre of the distribution at conditioning points in the
#' tails.
#'
#' @references
#' Gnasso, A. (2026). *Semiparametric Inference for Conditional Shapley Feature Importance*.
#' arXiv preprint, doi:10.48550/arXiv.2609.10313
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 300; p <- 4
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("X", seq_len(p))
#' f_lin <- function(Xm) 2 * Xm[, 1] - Xm[, 2]
#' result <- sage_cv(f_lin, X, K = 3, M = 20, B = 10,
#'                   family_set = "parametric", seed = 1)
#' print(result)
#' }
#'
#' @seealso [shapCopula::sage_cv_gcop()] for the faster Gaussian-copula variant;
#'   [shapCopula::sage_cv_marginal()] for the interventional baseline.
#'
#' @export
sage_cv <- function(f_predict, X, K = 3L, M = 30L, B = 20L,
                    alpha = 0.05, family_set = "nonparametric",
                    cores = 1L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Bh <- max(2L, as.integer(B %/% 2))

  fold_fn <- function(train_idx, test_idx, X) {
    X_tr     <- X[train_idx, , drop = FALSE]
    X_test   <- X[test_idx, , drop = FALSE]
    f_X_test <- f_predict(X_test)
    # The coalitions visited along a Shapley permutation are its prefixes, which
    # are the suffixes of the reversed order. Fitting one D-vine per permutation
    # with that order therefore makes every coalition exactly conditionable via
    # the Rosenblatt transform, at the cost of M vine fits per fold rather than
    # one fit and an approximate sampler.
    factory <- function(perm) {
      vc_list <- fit_copula(X_tr, family_set = family_set, cores = cores,
                            order = rev(perm))
      function(S, x_S_mat, B_per_row)
        vine_cond_sampler_exact(vc_list, S, x_S_mat, B_per_row)
    }
    .debiased_psi(f_predict, X_test, f_X_test, factory, M, Bh)
  }

  res <- .sage_crossfit(f_predict, X, K, alpha, fold_fn, M = M)
  attr(res, "sampler") <- "vine copula (exact Rosenblatt)"
  res
}


# --- Exported: Marginal SAGE (interventional baseline) -----------------------

#' Marginal (Interventional) SAGE with Cross-Fitting
#'
#' Computes global feature importance using the **marginal** (interventional)
#' value function \eqn{v_\text{marg}(S; x) = E_{X_{-S}}[f(x_S, X_{-S})]},
#' where \eqn{X_{-S}} is drawn independently from its marginal empirical
#' distribution (bootstrap rows of training data), ignoring the conditional
#' dependence given \eqn{X_S}. This corresponds to the interventional SHAP
#' framework (Janzing et al. 2020).
#'
#' @inheritParams sage_cv_gcop
#'
#' @return A `data.frame` with the same columns as [shapCopula::sage_cv_gcop()].
#'
#' @details
#' This estimator serves as a competitor / baseline for [shapCopula::sage_cv_gcop()].
#' When features are independent, both estimators coincide. Under dependence,
#' the marginal value function breaks the joint distribution and can
#' misattribute importance to spuriously correlated features.
#'
#' The marginal sampler is computationally inexpensive (no copula fitting)
#' and scales to arbitrary p.
#'
#' @references
#' Janzing, D., Minorics, L., & Bloebaum, P. (2020). Feature relevance
#' quantification in explainable AI: A causal problem.
#' *AISTATS*, 2907-2916.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 400; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("X", seq_len(p))
#' f_lin <- function(Xm) 2 * Xm[, 1] - Xm[, 2]
#' result <- sage_cv_marginal(f_lin, X, K = 3, M = 20, B = 10, seed = 1)
#' print(result)
#' }
#'
#' @seealso [shapCopula::sage_cv_gcop()] for the conditional (copula-based) estimator.
#'
#' @export
sage_cv_marginal <- function(f_predict, X, K = 3L, M = 30L, B = 20L,
                             alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Bh <- max(2L, as.integer(B %/% 2))

  fold_fn <- function(train_idx, test_idx, X) {
    X_tr     <- X[train_idx, , drop = FALSE]
    X_test   <- X[test_idx, , drop = FALSE]
    f_X_test <- f_predict(X_test)
    sampler  <- function(S, x_S_mat, B_per_row)
      marg_cond_sampler(X_tr, S, x_S_mat, B_per_row)
    .debiased_psi(f_predict, X_test, f_X_test, function(perm) sampler, M, Bh)
  }

  res <- .sage_crossfit(f_predict, X, K, alpha, fold_fn, M = M)
  attr(res, "sampler") <- "marginal (interventional)"
  res
}
