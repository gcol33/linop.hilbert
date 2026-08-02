## Closed-form ground truth, and the explicit matrix the matrix-free block must
## agree with. A closed form is code: everything here is checked against its own
## definition in test-closed-forms.R before any solver runs on it.

## V = v delta_0 has exactly one eigenvalue outside [-2, 2], with eigenvector
## u_j = rho^|j|. Substituting into (H u)_j = u_{j-1} + u_{j+1} + V_j u_j:
##   j = 0:  2 rho + v      = lambda            (using u_{-1} = u_1 = rho)
##   j > 0:  rho^{j-1} + rho^{j+1} = lambda rho^j
## The second gives rho + 1/rho = lambda, the first fixes the branch.
single_well_eigenvalue <- function(v) sign(v) * sqrt(v^2 + 4)
single_well_rho <- function(v) (-abs(v) + sqrt(v^2 + 4)) / 2
single_well_vector <- function(v, n) single_well_rho(v)^abs(seq.int(-n, n))

## H_n as a stored matrix, built from the definition rather than from the package.
explicit_block <- function(op, n) {
  N <- 2L * n + 1L
  M <- matrix(0, N, N)
  M[cbind(op$at + n + 1L, op$at + n + 1L)] <- op$V
  for (i in seq_len(N - 1L)) {
    M[i, i + 1L] <- 1
    M[i + 1L, i] <- 1
  }
  M
}

## The top of the spectrum of the block, by dense decomposition. Independent of
## linop::eigs(), so a test comparing the two is comparing two implementations.
dense_top <- function(op, n) {
  e <- eigen(explicit_block(op, n), symmetric = TRUE)
  i <- which.max(abs(e$values))
  list(value = e$values[i], vector = e$vectors[, i])
}
