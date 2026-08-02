## The operator this package certifies:
##
##     (H u)_j = u_{j-1} + u_{j+1} + V_j u_j        on l^2(Z),
##
## with V real and supported on finitely many sites. Two properties of that class
## are what make a closed-form certificate possible, and both are lost as soon as
## the support is infinite:
##
## 1. V finitely supported implies V compact, so by Weyl's theorem H and the free
##    operator have the same essential spectrum. The free operator is the Laurent
##    operator with symbol 2 cos(theta), so
##
##        sigma_ess(H) = [-2, 2].
##
##    The band edge is known a priori from the limiting coefficients rather than
##    estimated, and every bound below divides by a distance to it.
##
## 2. Outside the support the eigenvalue equation is the free recurrence
##    u_{j-1} + u_{j+1} = lambda u_j, whose characteristic roots are reciprocal.
##    For |lambda| > 2 exactly one lies in the open unit disc, so an l^2
##    eigenfunction is EXACTLY A rho^j out there, not asymptotically. That is what
##    turns the truncation error into a closed form rather than an existence
##    statement.
##
## No class is named `hilbert`. linop() is an S3 generic whose adapter convention
## is linop.<class>(), so a class of that name would collide with this package's
## own name.

#' A discrete Schrodinger operator on l^2(Z)
#'
#' The Jacobi operator `(H u)[j] = u[j-1] + u[j+1] + V[j] u[j]`, with a real
#' potential `V` supported on finitely many integer sites.
#'
#' The operator is self-adjoint by construction: the off-diagonals are `1` and the
#' diagonal `V` is real. Its essential spectrum is `[-2, 2]` whatever `V` is,
#' because a finitely supported potential is a compact perturbation of the free
#' operator, so anything outside that interval is an isolated eigenvalue of finite
#' multiplicity. Those are the values [eigenpairs()] computes and certifies.
#'
#' @param V The potential, a finite real vector. `V = 0` is allowed and describes
#'   the free operator, which has no eigenvalues at all; [eigenpairs()] reports
#'   that rather than inventing one.
#' @param at The integer sites `V` sits on, in the same order. Defaults to
#'   centring the support on the origin.
#' @return An object of class `finite_section`.
#' @seealso [discretise()] for the finite block, [eigenpairs()] for the certified
#'   eigenvalues.
#' @examples
#' ## a single well of depth 1 at the origin, whose eigenvalue is sqrt(5)
#' finite_section(1)
#'
#' ## a three-site potential
#' finite_section(c(1.5, -0.5, 2), at = -1:1)
#' @export
finite_section <- function(V, at = NULL) {
  ## Ahead of the numeric test, because a complex V is the interesting refusal:
  ## every bound this package issues rests on H being self-adjoint, and a complex
  ## diagonal breaks that with nothing downstream noticing.
  if (is.complex(V)) stopf("V must be real; a complex diagonal is not self-adjoint")
  if (!is.numeric(V) || !length(V) || !all(is.finite(V))) {
    stopf("V must be a non-empty finite numeric vector")
  }
  at <- at %||% (seq_along(V) - (length(V) + 1L) %/% 2L)
  if (!is_whole(at) || length(at) != length(V)) {
    stopf("at must be %d integer sites, one per entry of V", length(V))
  }
  at <- as.integer(at)
  if (anyDuplicated(at)) stopf("at repeats the site %d", at[anyDuplicated(at)])
  ord <- order(at)

  structure(list(V = as.numeric(V)[ord], at = at[ord],
                 support_radius = max(abs(at)),
                 band_edge = 2,
                 ## ||H|| <= ||free|| + ||V||, and both are exact for this class:
                 ## the free operator has norm 2 and V is a diagonal.
                 norm_bound = 2 + max(abs(V))),
            class = "finite_section")
}

#' @export
print.finite_section <- function(x, ...) {
  cat(sprintf("<finite_section> discrete Schrodinger operator on l^2(Z)\n"))
  cat(sprintf("  V supported on %d site%s, j in [%d, %d]\n",
              length(x$V), if (length(x$V) == 1L) "" else "s",
              min(x$at), max(x$at)))
  cat(sprintf("  V = %s\n", paste(format(x$V, digits = 5), collapse = "  ")))
  cat(sprintf("  essential spectrum [%g, %g], ||H|| <= %g\n",
              -x$band_edge, x$band_edge, x$norm_bound))
  invisible(x)
}

## ---------------------------------------------------------------- closed forms

#' The decay rate of an eigenfunction
#'
#' Outside the support of `V` an eigenfunction of [finite_section()] at `lambda`
#' is exactly `A rho^j`, where `rho` is the root of `z^2 - lambda z + 1` inside
#' the unit disc. This returns that `rho`.
#'
#' It is defined only outside the essential spectrum. For `|lambda| <= 2` both
#' roots lie on the unit circle, nothing decays, and there is no eigenfunction to
#' describe; `NA` is returned rather than a number that would read as one.
#'
#' @param lambda A numeric vector of values.
#' @return `|lambda|/2 - sqrt((lambda/2)^2 - 1)`, in `(0, 1)`, or `NA` where
#'   `|lambda| <= 2`.
#' @examples
#' ## the single-well eigenvalue sqrt(5) decays at the golden-ratio conjugate
#' decay_rate(sqrt(5))
#' @export
decay_rate <- function(lambda) {
  a <- abs(lambda) / 2
  ifelse(a > 1, a - sqrt(pmax(a^2 - 1, 0)), NA_real_)
}
