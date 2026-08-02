## The link back from the finite block to the operator it approximates.
##
## linop carries an opaque envelope and never inspects it; the four generics below
## dispatch on the payload's class, which is this package's. The payload is the
## record of a discretisation -- the original operator and the width it was cut at
## -- rather than the operator alone, because everything a caller of the finite
## object wants to ask needs both.

new_discretisation_record <- function(operator, n) {
  structure(list(operator = operator, n = n),
            class = "finite_section_provenance")
}

#' Provenance methods for a finite section
#'
#' Methods on [linop::provenance_lift()], [linop::provenance_refine()],
#' [linop::provenance_original_residual()] and [linop::provenance_summary()] for
#' the envelope [discretise()] attaches. A caller holding only the finite `linop`
#' reaches all four through [linop::provenance()].
#'
#' @param p A provenance envelope obtained from [linop::provenance()].
#' @param x A vector or matrix on the finite block.
#' @param n_new A new half-width.
#' @param result A matrix of eigenvectors on the finite block, one per column, or
#'   a single vector.
#' @param ... Unused.
#' @name finite_section_provenance
#' @examples
#' A <- discretise(finite_section(1), 8)
#' p <- linop::provenance(A)
#' linop::provenance_summary(p)
#' head(linop::provenance_lift(p, rep(1, 17)))
NULL

#' @rdname finite_section_provenance
#' @exportS3Method linop::provenance_summary
provenance_summary.finite_section_provenance <- function(p, ...) {
  op <- p$payload$operator
  sprintf("finite section at n = %d of a discrete Schrodinger operator, V on %d site%s in [%d, %d]",
          p$payload$n, length(op$V), if (length(op$V) == 1L) "" else "s",
          min(op$at), max(op$at))
}

#' @rdname finite_section_provenance
#' @exportS3Method linop::provenance_refine
provenance_refine.finite_section_provenance <- function(p, n_new, ...) {
  discretise(p$payload$operator, n_new)
}

#' @rdname finite_section_provenance
#' @exportS3Method linop::provenance_lift
provenance_lift.finite_section_provenance <- function(p, x, ...) {
  n <- p$payload$n
  x <- as.matrix(x)
  if (nrow(x) != 2L * n + 1L) {
    stopf("x has %d rows; the block at n = %d has %d sites",
          nrow(x), n, 2L * n + 1L)
  }
  ## The extension by zero, which is the embedding l^2(-n..n) -> l^2(Z). Every
  ## site outside the block is exactly zero, so the index is the whole content of
  ## the lift and no values are invented.
  data.frame(j = sites(n), x, check.names = FALSE)
}

#' @rdname finite_section_provenance
#' @exportS3Method linop::provenance_original_residual
provenance_original_residual.finite_section_provenance <- function(p, result, ...) {
  U <- as.matrix(result)
  n <- p$payload$n
  if (nrow(U) != 2L * n + 1L) {
    stopf("result has %d rows; the block at n = %d has %d sites",
          nrow(U), n, 2L * n + 1L)
  }
  ## Extend u by zero and apply H. Every site with |j| <= n cancels, because the
  ## couplings P_n removed multiply entries that are zero; every site with
  ## |j| > n+1 is zero term by term. What is left sits on exactly two sites,
  ##
  ##     r_{n+1} = u_n,    r_{-(n+1)} = u_{-n},
  ##
  ## so the residual against the ORIGINAL operator is two entries of the computed
  ## vector. No apply, no norm estimate, and nothing that grows with n.
  N <- nrow(U)
  sqrt(U[1L, ]^2 + U[N, ]^2) / sqrt(colSums(U^2))
}
