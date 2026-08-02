test_that("the matrix-free block is the explicit tridiagonal, entry for entry", {
  for (op in list(finite_section(1),
                  finite_section(c(1.5, -0.5, 2), at = -1:1),
                  finite_section(0),
                  finite_section(c(-3, 4), at = c(-2, 3)))) {
    n <- max(op$support_radius + 1L, 6L)
    A <- discretise(op, n)
    expect_equal(as.matrix(A), explicit_block(op, n),
                 info = paste(format(op$V), collapse = " "))
    expect_equal(dim(A), c(2L * n + 1L, 2L * n + 1L))
  }
})

test_that("the block satisfies linop's operator contract", {
  ## The claim the whole package rests on is that this object behaves like the
  ## self-adjoint operator it says it is. verify() is what checks that, and it is
  ## checking a callback here, not a stored matrix.
  A <- discretise(finite_section(c(1.5, -0.5, 2)), 12)
  cert <- linop::verify(A)
  expect_equal(cert$overall, "pass")
  expect_true(isTRUE(linop::cap(A, "hermitian")$value))
  expect_true(isTRUE(linop::cap(A, "symmetric")$value))
})

test_that("a block applies to several columns at once, and agrees column by column", {
  A <- discretise(finite_section(c(1, -1)), 9)
  X <- matrix(rnorm(nrow(as.matrix(A)) * 3), ncol = 3)
  Y <- A %*% X
  for (j in 1:3) {
    expect_equal(as.numeric(Y[, j]), as.numeric(A %*% X[, j, drop = FALSE]))
  }
})

test_that("a block narrower than the support is refused, with the reason", {
  op <- finite_section(c(1.5, -0.5, 2), at = -1:1)
  expect_error(discretise(op, 1), "does not exceed the support radius")
  expect_error(discretise(op, 0), "does not exceed the support radius")
  expect_error(discretise(op, -3), "does not exceed the support radius")
  ## the smallest block that is certifiable is the first one past the support
  expect_silent(discretise(op, 2))
})

test_that("discretise validates n", {
  op <- finite_section(1)
  expect_error(discretise(op, 3.5), "single integer")
  expect_error(discretise(op, c(4, 5)), "single integer")
  expect_error(discretise(op, "4"), "single integer")
})

test_that("sites indexes the block from -n to n", {
  expect_equal(sites(3), -3:3)
  expect_equal(length(sites(20)), 41L)
  expect_error(sites(2.5), "single integer")
})
