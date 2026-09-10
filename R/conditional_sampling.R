# Conditional samplers — internal workhorses for the SAGE estimators.
# Each sampler has a uniform interface:
#
#   sampler(S, x_S_matrix, B_per_row) -> matrix of shape (N * B_per_row) x p
#
# where N = nrow(x_S_matrix). For S = integer(0), x_S_matrix is ignored
# (pass matrix(0, 1, 0)) and the function returns N * B_per_row marginal draws.
# Rows are ordered as: first B_per_row rows belong to observation 1, next
# B_per_row to observation 2, and so on.


# --- Gaussian copula sampler -------------------------------------------------

# Sample X_{-S} | X_S = x_S under the Gaussian copula:
#   x_S -> u_S = F_S(x_S) -> z_S = Phi^{-1}(u_S)
#   z_{-S} | z_S ~ N(mu_cond, Sigma_cond)  [analytic]
#   x_{-S} = F_{-S}^{-1}(Phi(z_{-S}))      [quantile back-transform]
gcop_cond_sampler <- function(gcop, S, x_S_matrix, B) {
  X_train <- gcop$X_train
  n_tr    <- nrow(X_train)
  p       <- ncol(X_train)
  Sigma   <- gcop$Sigma
  N       <- if (length(S) == 0) 1L else nrow(x_S_matrix)

  if (length(S) == 0) {
    L <- chol(Sigma)
    Z <- matrix(stats::rnorm(N * B * p), N * B, p) %*% L
    U <- stats::pnorm(Z)
    X <- matrix(NA_real_, N * B, p)
    colnames(X) <- colnames(X_train)
    for (j in seq_len(p))
      X[, j] <- stats::quantile(gcop$sorted_cols[[j]], probs = U[, j],
                                names = FALSE, type = 1, na.rm = TRUE)
    return(X)
  }

  if (length(S) == p)
    return(as.matrix(x_S_matrix[rep(seq_len(N), each = B), , drop = FALSE]))

  not_S        <- setdiff(seq_len(p), S)
  Sigma_SS_inv <- solve(Sigma[S, S, drop = FALSE] + diag(1e-9, length(S)))
  Sigma_NS     <- Sigma[not_S, S, drop = FALSE]
  cond_cov     <- Sigma[not_S, not_S, drop = FALSE] -
                  Sigma_NS %*% Sigma_SS_inv %*% t(Sigma_NS)
  cond_cov     <- (cond_cov + t(cond_cov)) / 2
  L_cond <- tryCatch(
    chol(cond_cov),
    error = function(e) chol(cond_cov + diag(1e-7, nrow(cond_cov)))
  )

  out <- matrix(NA_real_, N * B, p)
  colnames(out) <- colnames(X_train)
  for (i in seq_len(N)) {
    xi_S <- as.numeric(x_S_matrix[i, ])
    u_S  <- numeric(length(S))
    for (k in seq_along(S)) {
      e      <- gcop$ecdfs[[S[k]]](xi_S[k])
      u_S[k] <- pmin(pmax(e, 1 / (n_tr + 1)), n_tr / (n_tr + 1))
    }
    z_S       <- stats::qnorm(u_S)
    cond_mean <- as.numeric(Sigma_NS %*% Sigma_SS_inv %*% z_S)
    Z_NS      <- matrix(stats::rnorm(B * length(not_S)), B, length(not_S)) %*% L_cond
    Z_NS      <- sweep(Z_NS, 2, cond_mean, "+")
    U_NS      <- stats::pnorm(Z_NS)
    X_NS      <- matrix(NA_real_, B, length(not_S))
    for (k in seq_along(not_S))
      X_NS[, k] <- stats::quantile(gcop$sorted_cols[[not_S[k]]], probs = U_NS[, k],
                                   names = FALSE, type = 1, na.rm = TRUE)
    rows        <- ((i - 1L) * B + 1L):(i * B)
    out[rows, S]     <- matrix(xi_S, B, length(S), byrow = TRUE)
    out[rows, not_S] <- X_NS
  }
  out
}


# --- Marginal (interventional) sampler --------------------------------------

# Sample X_{-S} | do(x_S): each out-of-coalition coordinate is drawn
# independently from its marginal empirical distribution (bootstrap rows of
# X_train), and x_S is held fixed. This implements the interventional / TreeSHAP
# value function (Janzing et al. 2020).
marg_cond_sampler <- function(X_train, S, x_S_matrix, B) {
  n_tr <- nrow(X_train)
  p    <- ncol(X_train)
  N    <- if (length(S) == 0) 1L else nrow(x_S_matrix)
  not_S <- setdiff(seq_len(p), S)

  out <- matrix(NA_real_, N * B, p)
  colnames(out) <- colnames(X_train)

  if (length(S) == 0) {
    for (k in seq_len(p))
      out[, k] <- X_train[sample.int(n_tr, N * B, replace = TRUE), k]
    return(out)
  }
  if (length(S) == p)
    return(as.matrix(x_S_matrix[rep(seq_len(N), each = B), , drop = FALSE]))

  for (i in seq_len(N)) {
    rows <- ((i - 1L) * B + 1L):(i * B)
    out[rows, S] <- matrix(as.numeric(x_S_matrix[i, ]), B, length(S), byrow = TRUE)
    for (k in not_S)
      out[rows, k] <- X_train[sample.int(n_tr, B, replace = TRUE), k]
  }
  out
}


# --- Vine copula conditional sampler (exact, via Rosenblatt) -----------------

# Sample X_{-S} | X_S = x_S exactly, using the Rosenblatt transform of a vine
# whose structure is a D-vine with S occupying the *last* positions of the
# variable order.
#
# For such a vine the Rosenblatt map w = R(u) is sequential from the end of the
# order, so the components w_S depend on u_S alone. Fixing those components and
# redrawing the remaining ones as independent uniforms therefore yields exact
# draws from the conditional law: no kernel smoothing, no importance weights.
#
# Callers must guarantee that S is a suffix of the order used to fit `vc`;
# `vine_perm_sampler` below arranges this by fitting one D-vine per Shapley
# permutation, whose nested coalitions are exactly the suffixes of that order.
vine_cond_sampler_exact <- function(vc_list, S, x_S_matrix, B) {
  p <- ncol(vc_list$X_train)
  N <- if (length(S) == 0) 1L else nrow(x_S_matrix)

  if (length(S) == 0) {
    U <- rvinecopulib::rvinecop(N * B, vc_list$vc)
    out <- from_pseudo_obs(U, vc_list$sorted_cols)
    colnames(out) <- colnames(vc_list$X_train)
    return(out)
  }
  if (length(S) == p)
    return(as.matrix(x_S_matrix[rep(seq_len(N), each = B), , drop = FALSE]))

  # x_S -> u_S on the copula scale (clipped away from the boundary).
  U_S <- to_pseudo_obs_cols(x_S_matrix, vc_list$ecdfs[S], vc_list$n_train)

  # w_S = Rosenblatt(u)_S depends only on u_S, so the remaining coordinates of
  # the input are irrelevant and are filled with 0.5.
  U_probe <- matrix(0.5, N, p)
  U_probe[, S] <- U_S
  W_S <- rvinecopulib::rosenblatt(U_probe, vc_list$vc)[, S, drop = FALSE]

  # One batched inverse transform for all N * B draws.
  W <- matrix(stats::runif(N * B * p), N * B, p)
  W[, S] <- W_S[rep(seq_len(N), each = B), , drop = FALSE]
  U_new <- rvinecopulib::inverse_rosenblatt(W, vc_list$vc)

  out <- from_pseudo_obs(U_new, vc_list$sorted_cols)
  # Pin the conditioning coordinates to the exact requested values: the
  # round-trip through the empirical quantile function is a step map and would
  # otherwise return the nearest training value.
  out[, S] <- as.matrix(x_S_matrix)[rep(seq_len(N), each = B), , drop = FALSE]
  colnames(out) <- colnames(vc_list$X_train)
  out
}
