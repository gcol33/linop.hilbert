## The envelope is obtained from discretise(), which is the only way a caller ever
## gets one. Hand-building it would test the generics and not the coupling, and
## would pass while the caller's path was broken.

test_that("all four generics route from an envelope a caller actually holds", {
  op <- finite_section(c(1.5, -0.5, 2), at = -1:1)
  A <- discretise(op, 10)
  p <- linop::provenance(A)

  expect_equal(p$provider, "linop.hilbert")
  expect_match(linop::provenance_summary(p), "finite section at n = 10")
  expect_match(linop::provenance_summary(p), "3 sites in \\[-1, 1\\]")

  B <- linop::provenance_refine(p, 25)
  expect_equal(dim(B), c(51L, 51L))
  expect_equal(as.matrix(B), explicit_block(op, 25))
  ## the refined operator carries its own envelope, so refinement composes
  expect_match(linop::provenance_summary(linop::provenance(B)), "n = 25")

  lifted <- linop::provenance_lift(p, seq_len(21))
  expect_equal(lifted$j, -10:10)
  expect_equal(lifted[[2L]], seq_len(21))
})

test_that("the envelope survives the algebra, so a composed operator still knows", {
  A <- discretise(finite_section(1), 6)
  expect_match(linop::provenance_summary(linop::provenance(t(A))), "n = 6")
  expect_match(linop::provenance_summary(linop::provenance(2 * A)), "n = 6")
  expect_match(linop::provenance_summary(linop::provenance(t(A) %*% A)), "n = 6")
})

test_that("the residual against the original operator is exactly the two boundary entries", {
  op <- finite_section(1)
  n <- 12L
  A <- discretise(op, n)
  p <- linop::provenance(A)

  ## Take any vector, not an eigenvector: the identity being asserted is about the
  ## extension by zero and holds for every u.
  set.seed(4)
  u <- matrix(rnorm(2L * n + 1L), ncol = 1L)
  reported <- linop::provenance_original_residual(p, u)

  ## the same quantity computed on a much wider block, where (H - 0) u~ can be
  ## formed explicitly and H is the original operator to within a truncation that
  ## touches nothing
  wide <- 40L
  u_ext <- numeric(2L * wide + 1L)
  u_ext[(wide - n + 1L):(wide + n + 1L)] <- u
  Hu <- explicit_block(op, wide) %*% u_ext
  ## r = (H - theta) u~ with theta = 0 is H u~; the claim is that outside the block
  ## it is supported on exactly two sites
  outside <- abs(sites(wide)) > n
  expect_equal(sum(Hu[outside] != 0), 2L)
  expect_equal(sqrt(sum(Hu[outside]^2)) / sqrt(sum(u^2)), reported,
               tolerance = 1e-14)
})

test_that("the residual is per column for a block of vectors", {
  A <- discretise(finite_section(1), 8)
  p <- linop::provenance(A)
  U <- cbind(c(1, rep(0, 15), 1), c(rep(0, 16), 3))
  eta <- linop::provenance_original_residual(p, U)
  expect_length(eta, 2L)
  expect_equal(eta, c(sqrt(2) / sqrt(2), 3 / 3))
})

test_that("a vector of the wrong length is refused rather than recycled", {
  p <- linop::provenance(discretise(finite_section(1), 5))
  expect_error(linop::provenance_lift(p, 1:5), "the block at n = 5 has 11 sites")
  expect_error(linop::provenance_original_residual(p, matrix(1, 4, 1)),
               "the block at n = 5 has 11 sites")
})
