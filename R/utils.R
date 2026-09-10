# Internal utilities shared across the package.

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Transform pseudo-observations U in (0,1)^p back to the original scale
# via empirical quantile (type-1 step function, ties = first).
from_pseudo_obs <- function(U, sorted_cols) {
  U <- as.matrix(U)
  p <- ncol(U)
  X <- matrix(NA_real_, nrow(U), p)
  for (j in seq_len(p))
    X[, j] <- stats::quantile(sorted_cols[[j]], probs = U[, j],
                              names = FALSE, type = 1, na.rm = TRUE)
  X
}

# Transform X (original scale) to pseudo-observations u_j = F_j(x_j) in (0,1),
# clipped to [1/(N+1), N/(N+1)] to avoid boundary values.
to_pseudo_obs <- function(X, ecdfs) {
  X <- as.matrix(X)
  p <- ncol(X)
  U <- matrix(NA_real_, nrow(X), p)
  for (j in seq_len(p)) {
    u <- ecdfs[[j]](X[, j])
    N <- length(environment(ecdfs[[j]])$x)
    U[, j] <- pmin(pmax(u, 1 / (N + 1)), N / (N + 1))
  }
  U
}

# Silverman bandwidth for a d-dimensional kernel on n observations.
silverman_bw <- function(n, d) {
  (4 / (d + 2))^(1 / (d + 4)) * n^(-1 / (d + 4))
}

# Transform selected columns of X to pseudo-observations, clipped to
# [1/(n+1), n/(n+1)]. `ecdfs_sub` holds the ECDFs of exactly those columns.
to_pseudo_obs_cols <- function(X, ecdfs_sub, n_train) {
  X <- as.matrix(X)
  U <- matrix(NA_real_, nrow(X), ncol(X))
  for (k in seq_along(ecdfs_sub))
    U[, k] <- pmin(pmax(ecdfs_sub[[k]](X[, k]), 1 / (n_train + 1)),
                   n_train / (n_train + 1))
  U
}
