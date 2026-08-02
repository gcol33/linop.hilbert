## The verb. One front door: give it the operator and a width, get back the values
## that survived certification, the full table of what was computed, and the
## certificate that decided between them.

#' Certified eigenvalues of an operator on l^2(Z)
#'
#' Truncates `x` to the block `|j| <= n`, computes the extreme eigenpairs of that
#' block with [linop::eigs()], and certifies each one against the spectrum of the
#' original operator.
#'
#' Every isolated eigenvalue of [finite_section()] lies outside `[-2, 2]` and
#' everything inside that interval is essential spectrum, so the values worth
#' asking for are the ones of largest magnitude. That is the default and the only
#' ordering that finds them all.
#'
#' A value is certified when it lies outside the essential spectrum. A value
#' inside it is a discretisation of continuous spectrum: the original operator has
#' no eigenvalue anywhere near it, and no residual makes it one. Those values are
#' returned, in `ritz`, and are not put in `values`, where they would read as
#' eigenvalues of `x`.
#'
#' @param x A [finite_section()].
#' @param n Half-width of the block. Larger `n` costs `O(n)` per apply and buys a
#'   bound that falls like [decay_rate()] to the power `n`.
#' @param k How many values to compute.
#' @param which Which end of the block's spectrum, passed to [linop::eigs()].
#'   `"largest"` is largest magnitude.
#' @param tol The tolerance the certificate is read against. It is the target for
#'   the distance to the spectrum of the original operator, not for the inner
#'   eigensolve, which is run tighter.
#' @param ncv Width of the Krylov subspace the block eigensolve works in. The
#'   default is wider than [linop::eigs()]'s own, because an eigenvalue of this
#'   class sits just outside a band of essential spectrum and the block's
#'   remaining values cluster against it. Measured on `finite_section(0.3)` at
#'   `n = 80`, `ncv = 21` converges nothing and lands 7.3e-04 from the truth,
#'   where `ncv = 40` reaches 2.7e-12; wider than 40 changed no digit in any case
#'   tried.
#' @param ... Passed to [linop::eigs()].
#' @return An object of class `finite_section_eigen` with elements
#'   \describe{
#'     \item{`values`}{The certified eigenvalues, in the order computed. Empty
#'       when none were.}
#'     \item{`ritz`}{Every value computed, with its residual against the original
#'       operator, its bound, its distance from the band edge, its bracket, and
#'       whether it was certified.}
#'     \item{`vectors`}{The computed vectors on the block, one per column, indexed
#'       by [sites()].}
#'     \item{`certificate`}{A [linop::build_certificate()] object.}
#'   }
#' @seealso [finite_section()], [discretise()].
#' @examples
#' ## a single well of depth 1, whose eigenvalue is sqrt(5) in closed form
#' H <- finite_section(1)
#' fit <- eigenpairs(H, n = 20)
#' fit
#' fit$values - sqrt(5)
#'
#' ## the free operator has no eigenvalues, and the certificate says so
#' eigenpairs(finite_section(0), n = 20)
#' @export
eigenpairs <- function(x, n, k = 1L, which = "largest", tol = 1e-8,
                       ncv = NULL, ...) {
  if (!inherits(x, "finite_section")) stopf("eigenpairs() expects a finite_section")
  k <- as.integer(k)
  if (is.na(k) || k < 1L) stopf("k must be a positive integer")
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) {
    stopf("tol must be a single positive number")
  }
  A <- discretise(x, n)
  n <- as.integer(n)
  ncv <- as.integer(ncv %||% min(2L * n + 1L, max(40L, 4L * k + 20L)))

  ## The inner solve is run at the arithmetic floor rather than at tol. Its error
  ## enters the reported bound through the Rayleigh quotient, and there is nothing
  ## to gain by letting it be the binding term: the truncation is what tol is
  ## about, and it is the part that costs n to improve.
  inner <- linop::eigs(A, k = k, which = which, ncv = ncv,
                       tol = min(tol, .Machine$double.eps^0.75), ...)

  U <- as.matrix(inner$vectors)
  theta <- as.numeric(inner$values)
  ## linop::eigs() reports the Rayleigh quotient rather than the Ritz value, which
  ## is what Kato-Temple needs and what minimises the residual over the value.
  eta <- linop::provenance_original_residual(linop::provenance(A), U)

  cert <- finite_section_certificate(x, n, theta, U, eta, tol, inner, ncv)
  v <- cert$values

  ritz <- data.frame(value = theta, eta = eta, bound = v$bound, gap = v$gap,
                     lower = v$bracket[, "lower"], upper = v$bracket[, "upper"],
                     certified = v$isolated)

  structure(list(values = theta[v$isolated], ritz = ritz, vectors = U,
                 sites = sites(n), n = n, k = k, operator = x,
                 certificate = cert, inner = inner),
            class = "finite_section_eigen")
}

#' @export
print.finite_section_eigen <- function(x, ...) {
  op <- x$operator
  cat(sprintf("<finite_section_eigen> n = %d, block of %d sites\n",
              x$n, 2L * x$n + 1L))
  nc <- sum(x$ritz$certified)
  cat(sprintf("  %d of %d value%s certified outside the essential spectrum [%g, %g]\n",
              nc, nrow(x$ritz), if (nrow(x$ritz) == 1L) "" else "s",
              -op$band_edge, op$band_edge))
  for (i in seq_len(nrow(x$ritz))) {
    row <- x$ritz[i, ]
    if (row$certified) {
      cat(sprintf("    %+.12g   dist to sigma(H) <= %.3e   in [%.12g, %.12g]\n",
                  row$value, row$bound, row$lower, row$upper))
    } else {
      cat(sprintf("    %+.12g   inside the band, not an eigenvalue of H\n", row$value))
    }
  }
  cat(sprintf("  certificate: %s\n", x$certificate$overall))
  invisible(x)
}

#' @export
summary.finite_section_eigen <- function(object, ...) {
  print(object)
  print(object$certificate)
  invisible(object)
}
