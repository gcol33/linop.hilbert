## The certificate this package issues, which is a shape of its own.
##
## linop's four shapes all certify an answer to a question posed on the object in
## hand: does this operator behave like a linear map, how close is x to solving
## A x = b, how close is (theta, x) to an eigenpair OF THIS OPERATOR. Here the
## operator in hand is not the one the question is about. theta is an eigenvalue
## of H_n and the claim is about sigma(H), so the rows have to carry the step
## between them, and that step is the whole content:
##
##   arithmetic floor       what no computation below can resolve
##   block eigensolve       whether the finite problem was solved at all
##   original residual      ||(H - theta) u~|| / ||u~||, exactly two entries
##   truncation bound       dist(theta, sigma(H)) <= eta + floor, by self-adjointness
##   isolation gap          q - a > 0, without which none of the above says
##                          "eigenvalue", because the value may be a discretisation
##                          of continuous spectrum
##
## The block row is a qualification and never a failure, and the reason is worth
## stating because it is what makes the whole scheme robust. The residual bound
##
##     dist(theta, sigma(H)) <= ||(H - theta) u|| / ||u||
##
## holds for EVERY theta and every nonzero u, not only for eigenpairs. So a block
## eigensolve that stalled still yields a true bound; it yields a worse one. What
## the row reports is that the number the caller is reading was limited by the
## inner solve rather than by the truncation, which is the case where raising n
## makes the answer worse instead of better.

## c in the arithmetic floor. Higham's bound on a computed residual counts terms in
## an inner product, and here the inner product is three-term whatever n is, since
## every row of H has three nonzeros. S0.6 measured the resulting plateau at about
## 2.2e-15 against ||H|| of 4 to 6, which puts c between 2 and 4. The same default
## linop uses, for the same measurement.
FLOOR_CONST <- 4

## Kato-Temple, with the band edge known rather than estimated. With q the Rayleigh
## quotient, eta the residual norm, and no spectrum of H in the open interval
## (a, q),
##
##     q - eta^2 / (b - q)  <=  lambda  <=  q + eta^2 / (q - a).
##
## For the eigenvalue above the band, a = 2 and there is nothing above it to play
## b, so b is +Inf and the lower end is q itself: a Rayleigh quotient never exceeds
## the largest eigenvalue. Below the band the picture is reflected. The half-width
## is eta^2/(q - a) against the linear bound's eta, which at eta = 4e-5 is four
## orders of magnitude.
kato_temple <- function(q, eta, band_edge) {
  gap <- abs(q) - band_edge
  width <- eta^2 / gap
  if (q > 0) c(lower = q, upper = q + width) else c(lower = q - width, upper = q)
}

## Self-adjointness is what turns a residual into a spectral statement, and for
## this class it is a property of the constructor rather than a declaration: the
## off-diagonals are 1 and finite_section() refuses a complex V. So the row that
## rests on it records construction, and evidence_satisfies() sees that underneath
## the bound the way linop's forward-error row sees whatever established hermitian.
self_adjoint_evidence <- function() {
  linop::evidence("construction", "identity", 1)
}

finite_section_certificate <- function(op, n, theta, U, eta, tol, inner, ncv) {
  r <- linop::cert_rows()
  eps <- .Machine$double.eps
  floor_abs <- FLOOR_CONST * op$norm_bound * eps
  a <- op$band_edge
  k <- length(theta)

  gap <- abs(theta) - a
  isolated <- gap > 0
  bound <- eta + floor_abs
  bracket <- lapply(seq_len(k), function(i) {
    if (isolated[i]) kato_temple(theta[i], eta[i], a) else c(lower = NA, upper = NA)
  })

  r$add("arithmetic floor", "pass",
        sprintf("||H|| <= 2 + max|V| = %g exactly; c = %g, eps = %.3g, floor %.3e",
                op$norm_bound, FLOOR_CONST, eps, floor_abs),
        source = "construction", guarantee = "identity")

  r$add("block eigensolve",
        if (inner$nconv >= inner$k) "pass" else "qualified",
        sprintf("%d of %d pairs of H_n converged in %d of at most %d iterations at ncv = %d; %s",
                inner$nconv, inner$k, inner$iterations, inner$maxit, ncv,
                if (inner$nconv >= inner$k) "the bound below is set by the truncation"
                else "the bound below is set by this rather than by the truncation, and a wider ncv is what moves it"),
        source = "computation", guarantee = "identity")

  ## Exact rather than estimated: u~ is u extended by zero, (H - theta) u~ vanishes
  ## at every site the block covers and at every site beyond n+1, and the two that
  ## survive are entries of u. Nothing here is measured by applying an operator.
  r$add("original residual",
        if (all(eta <= tol)) "pass" else "qualified",
        sprintf("worst ||(H - theta) u~|| / ||u~|| = %.3e against tol %.3e, from the sites j = +-%d",
                max(eta), tol, n + 1L),
        source = "computation", guarantee = "identity")

  ## H is self-adjoint, so the distance from any theta to its spectrum is at most
  ## the residual norm of any unit vector. That is a theorem, and it is the reason
  ## the row above is worth computing.
  r$add("truncation bound",
        if (all(bound <= tol)) "pass" else if (all(is.finite(bound))) "qualified" else "fail",
        sprintf("dist(theta, sigma(H)) <= %.3e for every value reported", max(bound)),
        evidence = linop::evidence("theorem", "deterministic_bound", 1,
                                   depends_on = list(self_adjoint_evidence())))

  ## The one row that can fail on a correct computation, and the reason the package
  ## refuses rather than reporting a number. sigma_ess(H) = [-2, 2] for every V, so
  ## a Ritz value inside it is a discretisation of continuous spectrum and there is
  ## no eigenvalue anywhere near it to be close to. The same inequality is what the
  ## Kato-Temple denominator needs, so isolation and refusal are one condition.
  r$add("isolation gap", if (all(isolated)) "pass" else "fail",
        if (all(isolated))
          sprintf("|theta| - %g >= %.3e for every value; sigma_ess(H) = [%g, %g] by Weyl",
                  a, min(gap), -a, a)
        else
          sprintf(paste0("%d of %d values lie inside the essential spectrum [%g, %g], ",
                         "where H has no eigenvalues; those are discretisations of ",
                         "continuous spectrum and are not certified"),
                  sum(!isolated), k, -a, a),
        source = "computation", guarantee = "identity")

  linop::build_certificate(
    r$collect(), subject = "finite section",
    values = list(eta = eta, floor = floor_abs, bound = bound, gap = gap,
                  bracket = do.call(rbind, bracket), isolated = isolated,
                  n = n, ncv = ncv, inner = inner$certificate),
    evidence = r$collect_evidence())
}
