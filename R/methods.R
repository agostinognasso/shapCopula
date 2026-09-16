# S3 methods for the object returned by the sage_cv_* estimators.

# The estimators return a plain data.frame carrying the Wald summary; tagging it
# with a class lets print() and plot() do something sensible without taking away
# the ability to treat it as an ordinary data.frame.

#' Print a SAGE Estimate
#'
#' @param x An object of class `sage_estimate`, as returned by
#'   [shapCopula::sage_cv_gcop()], [shapCopula::sage_cv()] or
#'   [shapCopula::sage_cv_marginal()].
#' @param digits Integer. Significant digits for the numeric columns.
#' @param ... Passed to [print.data.frame()].
#'
#' @return `x`, invisibly. Called for the side effect of printing.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 3), 200, 3)
#' colnames(X) <- paste0("X", 1:3)
#' res <- sage_cv_gcop(function(Xm) 2 * Xm[, 1], X, K = 2, M = 5, B = 6, seed = 1)
#' print(res)
#'
#' @export
print.sage_estimate <- function(x, digits = 3, ...) {
  cat("Conditional SAGE values with Wald inference\n")
  cat("  sampler: ", attr(x, "sampler") %||% "unknown",
      "   folds: ", attr(x, "K") %||% NA,
      "   permutations: ", attr(x, "M") %||% NA,
      "   level: ", 100 * (1 - (attr(x, "alpha") %||% 0.05)), "%\n\n", sep = "")
  # Defaults that the caller may override through `...`.  Passing them to
  # print.data.frame() alongside `...` would raise "formal argument matched by
  # multiple actual arguments" as soon as the caller supplied one of them.
  dots <- list(...)
  if (!"row.names" %in% names(dots)) dots$row.names <- FALSE
  do.call(print.data.frame,
          c(list(as.data.frame(x), digits = digits), dots))
  invisible(x)
}

#' Forest Plot of SAGE Estimates
#'
#' Draws one horizontal confidence interval per feature, ordered by estimated
#' importance, with a reference line at zero. Features whose interval excludes
#' zero are drawn in a solid style, the others dashed, so that the outcome of
#' the test is readable directly off the plot.
#'
#' @param x An object of class `sage_estimate`.
#' @param sort Logical. Order features by `psi_hat` (default `TRUE`).
#' @param xlab,main Axis label and title.
#' @param ... Passed to [graphics::plot.default()].
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(300 * 4), 300, 4)
#' colnames(X) <- paste0("X", 1:4)
#' f <- function(Xm) 2 * Xm[, 1] - Xm[, 2]
#' res <- sage_cv_gcop(f, X, K = 2, M = 10, B = 8, seed = 1)
#' plot(res)
#'
#' @export
plot.sage_estimate <- function(x, sort = TRUE, xlab = expression(hat(Psi)[j]),
                               main = "Conditional SAGE with Wald intervals",
                               ...) {
  d <- as.data.frame(x)
  if (sort) d <- d[order(d$psi_hat), , drop = FALSE]
  n <- nrow(d)
  sig <- (d$ci_lo > 0) | (d$ci_hi < 0)

  op <- graphics::par(mar = c(4.5, max(5, 0.6 * max(nchar(d$feature)) + 3), 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  xlim <- range(c(d$ci_lo, d$ci_hi, 0), finite = TRUE)
  xlim <- xlim + c(-1, 1) * 0.04 * diff(xlim)
  # Same reasoning as in print.sage_estimate(): these are defaults, not fixed
  # values, so they are dropped whenever `...` already carries them.
  dots     <- list(...)
  defaults <- list(xlim = xlim, ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "")
  defaults <- defaults[setdiff(names(defaults), names(dots))]
  do.call(graphics::plot,
          c(list(NA, xlab = xlab, main = main), defaults, dots))
  graphics::abline(v = 0, col = "grey60", lty = 3)
  graphics::axis(2, at = seq_len(n), labels = d$feature, las = 1, tick = FALSE)
  graphics::segments(d$ci_lo, seq_len(n), d$ci_hi, seq_len(n),
                     lwd = ifelse(sig, 2, 1), lty = ifelse(sig, 1, 2))
  graphics::points(d$psi_hat, seq_len(n), pch = ifelse(sig, 19, 21),
                   bg = "white", cex = 1.1)
  invisible(x)
}
