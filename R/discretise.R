## The finite section H_n = P_n H P_n on the sites |j| <= n, as a linop.
##
## Matrix-free rather than a stored tridiagonal, because the apply is three shifted
## adds and a stored matrix would be 2n+1 squared for no gain. The operator carries
## the original one as its provenance payload, which is what lets a caller of the
## finite object recover what it approximates.

#' Truncate an operator to a finite block
#'
#' @param x The operator.
#' @param n Half-width of the block.
#' @param ... Passed to methods.
#' @return A [linop::linop()].
#' @export
discretise <- function(x, n, ...) UseMethod("discretise")

#' @describeIn discretise Keep the sites `|j| <= n`, giving a block of dimension
#'   `2n + 1`. `n` must exceed the support radius of `V`, so that the truncated
#'   couplings all fall in the region where the free recurrence holds; that is the
#'   condition under which the truncation error is the closed form
#'   [eigenpairs()] reports and not an estimate.
#' @export
discretise.finite_section <- function(x, n, ...) {
  if (!is_whole(n) || length(n) != 1L) stopf("n must be a single integer")
  n <- as.integer(n)
  if (n <= x$support_radius) {
    stopf(paste0("n = %d does not exceed the support radius %d.\n",
                 "  The truncation bound holds because the sites the block drops are ones\n",
                 "  where V is already zero and the free recurrence holds exactly. Inside\n",
                 "  the support there is no closed form to certify against."),
          n, x$support_radius)
  }
  N <- 2L * n + 1L
  d <- numeric(N)
  d[x$at + n + 1L] <- x$V

  ## (H u)_j = u_{j-1} + u_{j+1} + V_j u_j, with the sites outside the block
  ## reading as zero, which is what P_n H P_n does.
  apply_block <- function(X) {
    X <- as.matrix(X)
    rbind(X[-1L, , drop = FALSE], 0) + rbind(0, X[-nrow(X), , drop = FALSE]) + d * X
  }

  A <- linop::linop(apply_block, adjoint = apply_block, dim = c(N, N),
                    properties = c("hermitian", "symmetric"))
  linop::set_provenance(A, PROVIDER, new_discretisation_record(x, n))
}

#' The sites a discretisation covers
#'
#' @param n Half-width of the block, as passed to [discretise()].
#' @return The integer vector `-n:n`, which indexes the rows of the block and the
#'   entries of an eigenvector returned by [eigenpairs()].
#' @examples
#' sites(3)
#' @export
sites <- function(n) {
  if (!is_whole(n) || length(n) != 1L) stopf("n must be a single integer")
  seq.int(-as.integer(n), as.integer(n))
}
