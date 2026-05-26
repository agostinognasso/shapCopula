#' Fit a Gaussian Copula with Empirical Marginals
#'
#' Estimates the joint distribution of `X` using a Gaussian copula: marginals
#' are modelled non-parametrically via the empirical CDF, and the copula
#' dependence structure is captured by the correlation matrix of the
#' normal-quantile-transformed pseudo-observations.
#'
#' @param X Numeric matrix or data frame of shape n-by-p.
#'
#' @return A list with components:
#' \describe{
#'   \item{`Sigma`}{p-by-p latent Gaussian correlation matrix.}
#'   \item{`ecdfs`}{List of p empirical CDF functions.}
#'   \item{`sorted_cols`}{List of p sorted column vectors (for quantile
#'     back-transformation).}
#'   \item{`X_train`}{The original training matrix (stored for use by the
#'     conditional sampler).}
#' }
#'
#' @details
#' Pseudo-observations are computed as `u_j = rank(X_j) / (n + 1)` to avoid
#' boundary values. The latent correlation matrix is the sample correlation of
#' `qnorm(U)`. Eigenvalue flooring at `1e-8` ensures positive definiteness.
#'
#' @seealso [shapCopula::sage_cv_gcop()] for the main estimation workflow.
#'
#' @references
#' Gnasso, A. (2026). *Inference for Conditional Shapley Values via Vine Copulas*.
#' Manuscript under review.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 4), 200, 4)
#' gc <- fit_gauss_copula(X)
#' dim(gc$Sigma)  # 4 x 4
#'
#' @export
fit_gauss_copula <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  ecdfs       <- lapply(seq_len(p), function(j) stats::ecdf(X[, j]))
  sorted_cols <- lapply(seq_len(p), function(j) sort(X[, j]))
  U <- matrix(NA_real_, n, p)
  for (j in seq_len(p))
    U[, j] <- rank(X[, j], ties.method = "average") / (n + 1)
  Z <- stats::qnorm(U)
  Sigma <- stats::cor(Z)
  e <- eigen(Sigma, symmetric = TRUE)
  if (any(e$values < 1e-8)) {
    e$values <- pmax(e$values, 1e-8)
    Sigma <- e$vectors %*% diag(e$values) %*% t(e$vectors)
    Sigma <- (Sigma + t(Sigma)) / 2
    diag(Sigma) <- 1
  }
  list(ecdfs = ecdfs, sorted_cols = sorted_cols, Sigma = Sigma, X_train = X)
}


#' Fit a Nonparametric Vine Copula with Empirical Marginals
#'
#' Estimates the joint distribution of `X` using a nonparametric vine copula
#' (transformation local-likelihood, TLL; Geenens 2014) fitted via
#' `rvinecopulib::vinecop`. Marginals are modelled via the empirical CDF.
#'
#' @param X Numeric matrix or data frame of shape n-by-p.
#' @param family_set Character. Copula family set passed to
#'   `rvinecopulib::vinecop`. Default `"nonparametric"` uses kernel TLL pairs.
#'   Use `"parametric"` or a specific family vector for faster fitting.
#' @param cores Integer. Number of cores for parallel pair-copula fitting
#'   (passed to `rvinecopulib::vinecop`).
#'
#' @return A list with components:
#' \describe{
#'   \item{`vc`}{A `vinecop` object from `rvinecopulib`.}
#'   \item{`ecdfs`}{List of p empirical CDF functions.}
#'   \item{`sorted_cols`}{List of p sorted column vectors.}
#'   \item{`X_train`}{The original training matrix.}
#' }
#'
#' @details
#' **Computational note.** Vine copulas scale as O(p^2) pair-copulas. For
#' nonparametric families this is feasible for p up to approximately 10-15 on
#' standard hardware. For larger p, use [shapCopula::fit_gauss_copula()] instead.
#'
#' @seealso [shapCopula::sage_cv()] for the vine-copula-based SAGE workflow.
#'   [shapCopula::fit_gauss_copula()] for a faster Gaussian-copula alternative.
#'
#' @references
#' Geenens, G., Charpentier, A., & Paindaveine, D. (2017).
#' Probit transformation for nonparametric kernel estimation of the copula
#' density. *Bernoulli*, 23(3), 1848-1873.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(150 * 3), 150, 3)
#' vc <- fit_copula(X, family_set = "parametric")
#' class(vc$vc)  # "vinecop"
#' }
#'
#' @export
fit_copula <- function(X, family_set = "nonparametric", cores = 1) {
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.matrix(X)
  p <- ncol(X)
  ecdfs       <- lapply(seq_len(p), function(j) stats::ecdf(X[, j]))
  sorted_cols <- lapply(seq_len(p), function(j) sort(X[, j]))
  U <- rvinecopulib::pseudo_obs(X, ties_method = "average")
  vc <- rvinecopulib::vinecop(
    U,
    family_set    = family_set,
    cores         = cores,
    nonpar_method = "constant"
  )
  list(vc = vc, ecdfs = ecdfs, sorted_cols = sorted_cols, X_train = X)
}
