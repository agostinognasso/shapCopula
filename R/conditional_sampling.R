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


# --- Vine copula batch sampler (importance-sampling) ------------------------

# Sample X_{-S} | X_S = x_S using importance-reweighted draws from the joint
# vine copula. A pool of N_pool draws is reweighted by a Gaussian kernel on
# the x_S coordinates; the resulting weighted sample is rounded to the
# conditioning value on columns S for exact enforcement.
vine_cond_sampler_batch <- function(vc_list, S, x_S_matrix, B) {
  p  <- ncol(vc_list$X_train)
  N  <- if (length(S) == 0) 1L else nrow(x_S_matrix)

  if (length(S) == 0) {
    U  <- rvinecopulib::rvinecop(N * B, vc_list$vc)
    return(from_pseudo_obs(U, vc_list$sorted_cols))
  }
  if (length(S) == p)
    return(as.matrix(x_S_matrix[rep(seq_len(N), each = B), , drop = FALSE]))

  n_tr      <- nrow(vc_list$X_train)
  pool_size <- max(2000L, N * B * 2L)
  U_pool    <- rvinecopulib::rvinecop(pool_size, vc_list$vc)
  X_pool    <- from_pseudo_obs(U_pool, vc_list$sorted_cols)

  d    <- length(S)
  sd_S <- apply(vc_list$X_train[, S, drop = FALSE], 2, stats::sd)
  bw   <- silverman_bw(n_tr, d) * sd_S  # per-coordinate bandwidths

  out    <- matrix(NA_real_, N * B, p)
  pool_S <- X_pool[, S, drop = FALSE]

  for (i in seq_len(N)) {
    z    <- sweep(pool_S, 2, as.numeric(x_S_matrix[i, ]), "-")
    z    <- sweep(z, 2, bw, "/")
    logw <- -0.5 * rowSums(z^2)
    w    <- exp(logw - max(logw))
    if (sum(w) < 1e-10) w <- rep(1, pool_size)
    idx    <- sample.int(pool_size, B, replace = TRUE, prob = w)
    chosen <- X_pool[idx, , drop = FALSE]
    chosen[, S] <- matrix(as.numeric(x_S_matrix[i, ]), B, d, byrow = TRUE)
    rows <- ((i - 1L) * B + 1L):(i * B)
    out[rows, ] <- chosen
  }
  out
}
