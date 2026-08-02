## A closed form is code, and a wrong one passes plausible tests. Each of these
## checks a formula against its own definition, not against another formula, and
## the eigenvector is checked by residual rather than by comparing sorted spectra:
## a sign error inside a cosine can leave an eigenvalue SET invariant while pairing
## every value with the wrong vector.

test_that("the single-well eigenpair satisfies the recurrence on l^2(Z)", {
  for (v in c(0.5, 1, 2, 4, -1, -3)) {
    lambda <- single_well_eigenvalue(v)
    rho <- single_well_rho(v)
    n <- 30L
    u <- single_well_vector(v, n)
    if (v < 0) u <- u * (-1)^abs(seq.int(-n, n))

    ## (H u)_j at every interior site, straight from the definition
    j <- seq.int(-n, n)
    Vj <- ifelse(j == 0L, v, 0)
    interior <- 2:(2L * n)
    Hu <- u[interior - 1L] + u[interior + 1L] + Vj[interior] * u[interior]
    expect_equal(Hu, lambda * u[interior], tolerance = 1e-12,
                 info = sprintf("v = %g", v))

    ## and the value is outside the essential spectrum, which is what makes it an
    ## eigenvalue rather than a discretisation of continuous spectrum
    expect_gt(abs(lambda), 2)
    expect_gt(rho, 0)
    expect_lt(rho, 1)
  }
})

test_that("decay_rate inverts the characteristic equation", {
  for (v in c(0.5, 1, 2, 4)) {
    lambda <- single_well_eigenvalue(v)
    rho <- decay_rate(lambda)
    ## rho is a root of z^2 - lambda z + 1
    expect_equal(rho^2 - lambda * rho + 1, 0, tolerance = 1e-13)
    ## and it is the one inside the disc, matching the closed form for this model
    expect_equal(rho, single_well_rho(v), tolerance = 1e-13)
  }
  ## the golden-ratio conjugate, which is the case worth recognising by eye
  expect_equal(decay_rate(sqrt(5)), (sqrt(5) - 1) / 2, tolerance = 1e-14)
})

test_that("decay_rate is NA inside the essential spectrum, where nothing decays", {
  expect_true(is.na(decay_rate(0)))
  expect_true(is.na(decay_rate(2)))
  expect_true(is.na(decay_rate(-1.999)))
  expect_equal(is.na(decay_rate(c(0, 3, -3, 1))), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("the explicit block is symmetric tridiagonal with the potential on it", {
  op <- finite_section(c(1.5, -0.5, 2), at = -1:1)
  M <- explicit_block(op, 5)
  expect_equal(M, t(M))
  expect_equal(diag(M), c(0, 0, 0, 0, 1.5, -0.5, 2, 0, 0, 0, 0))
  expect_true(all(M[row(M) - col(M) == 1] == 1))
  expect_true(all(M[abs(row(M) - col(M)) > 1] == 0))
})
