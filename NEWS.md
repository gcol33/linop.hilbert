# linop.hilbert 0.0.0.9000

* First unit: `finite_section()` for the discrete Schrodinger operator on
  `l^2(Z)`, `discretise()` for the finite block as a `linop`, and `eigenpairs()`
  for the certified isolated eigenvalues.
* The certificate carries four rows: the arithmetic floor, the residual against
  the original operator, the truncation bound, and the isolation gap. It is a
  shape of its own, since the operator computed on is not the one the claim is
  about.
* A value inside the essential spectrum is returned and is not certified. The
  free operator `V = 0` yields none, which is a certificate rather than an error.
