## Recovery against closed-form truth, run to convergence, across several
## potentials and widths. A shape test would pass on an operator that returned the
## wrong eigenvalue.

test_that("the single-well eigenvalue is recovered, and the bound holds at every n", {
  for (v in c(0.5, 1, 2, 4, -1, -3)) {
    truth <- single_well_eigenvalue(v)
    H <- finite_section(v)
    for (n in c(5L, 10L, 20L, 40L)) {
      fit <- eigenpairs(H, n = n)
      expect_length(fit$values, 1L)
      err <- abs(fit$values[1L] - truth)
      bound <- fit$ritz$bound[1L]
      expect_lte(err, bound)
      expect_true(fit$ritz$certified[1L],
                  info = sprintf("v = %g, n = %d", v, n))
    }
    ## and the bound is not vacuous: by n = 40 it has bought most of the digits
    expect_lt(abs(eigenpairs(H, n = 40L)$values[1L] - truth), 1e-8)
  }
})

test_that("the bound falls at the closed-form decay rate", {
  v <- 1
  rho <- single_well_rho(v)
  H <- finite_section(v)
  ns <- c(10L, 20L, 30L)
  eta <- vapply(ns, function(n) eigenpairs(H, n = n)$ritz$eta[1L], numeric(1))
  ## eta / rho^n is a constant of the potential, so consecutive ratios agree
  ratios <- eta / rho^ns
  expect_equal(ratios[1L], ratios[2L], tolerance = 1e-3)
  expect_equal(ratios[2L], ratios[3L], tolerance = 1e-3)
})

test_that("the Kato-Temple bracket contains the truth and beats the linear bound", {
  for (v in c(1, 2, -2)) {
    truth <- single_well_eigenvalue(v)
    for (n in c(5L, 10L, 20L)) {
      fit <- eigenpairs(finite_section(v), n = n)
      row <- fit$ritz[1L, ]
      expect_gte(truth, row$lower - 1e-14)
      expect_lte(truth, row$upper + 1e-14)
      ## the quadratic improvement is the reason the bracket is reported at all
      expect_lt(row$upper - row$lower, row$eta)
    }
  }
})

test_that("the arithmetic floor is load-bearing, and it can be shown", {
  ## Past some n the analytic residual has decayed below what the eigensolve can
  ## resolve. There eta stops measuring truncation, the true error can exceed it,
  ## and only the floor keeps the bound honest. Without it a fully converged answer
  ## certifies as wrong, which is the failure mode S0.6 found.
  ##
  ## Where that happens depends on how well the block was solved, so the n is not a
  ## constant: the sweep looks for the crossing rather than asserting it at a width
  ## that would drift with the eigensolver.
  truth <- single_well_eigenvalue(1)
  ns <- c(40L, 60L, 80L, 100L, 140L)
  bit <- FALSE
  for (n in ns) {
    fit <- eigenpairs(finite_section(1), n = n)
    err <- abs(fit$values[1L] - truth)
    eta <- fit$ritz$eta[1L]
    floor_v <- fit$certificate$values$floor

    ## the bound as reported holds at every width, which is the claim
    expect_lte(err, eta + floor_v, label = sprintf("n = %d", n))
    expect_equal(fit$ritz$bound[1L], eta + floor_v)
    if (err > eta) bit <- TRUE
  }
  ## and somewhere in that sweep the floor was the only thing making it hold
  expect_true(bit)
})

test_that("eta stops measuring truncation once it reaches the eigensolve's noise", {
  ## The other side of the same fact, and the reason raising n is not always an
  ## improvement: past the crossing eta is set by how well H_n was solved, so it
  ## stops falling like the decay rate and moves around at the 1e-16 level.
  truth <- single_well_eigenvalue(1)
  rho <- single_well_rho(1)
  ns <- c(20L, 40L, 140L)
  eta <- vapply(ns, function(n) eigenpairs(finite_section(1), n = n)$ritz$eta[1L],
                numeric(1))
  ratio <- eta / rho^ns

  ## while truncation is what eta measures, eta / rho^n is a constant of the
  ## potential, so it falls by four orders from n = 20 to n = 40
  expect_lt(eta[2L], eta[1L] * 1e-3)
  expect_equal(ratio[1L], ratio[2L], tolerance = 1e-3)

  ## by n = 140 the analytic residual is around 1e-29 and eta is nowhere near it:
  ## it has stopped following the decay law entirely
  expect_gt(ratio[3L], ratio[2L] * 1e6)
  ## and it is still small in absolute terms, so nothing is wrong with the answer
  expect_lt(eta[3L], 1e-11)
  expect_lt(abs(eigenpairs(finite_section(1), n = 140L)$values[1L] - truth), 1e-14)
})

test_that("a potential with no closed form is recovered against a dense reference", {
  op <- finite_section(c(1.5, -0.5, 2), at = -1:1)
  ref <- dense_top(op, 400L)$value
  expect_gt(abs(ref), 2)
  for (n in c(5L, 10L, 20L, 40L)) {
    fit <- eigenpairs(op, n = n)
    expect_lte(abs(fit$values[1L] - ref), fit$ritz$bound[1L])
  }
})

test_that("the block eigenvalue agrees with a dense decomposition of the same block", {
  ## Two implementations of the same finite problem: linop::eigs() on a callback,
  ## and eigen() on the stored matrix.
  op <- finite_section(c(1.5, -0.5, 2), at = -1:1)
  for (n in c(6L, 15L)) {
    fit <- eigenpairs(op, n = n)
    expect_equal(fit$ritz$value[1L], dense_top(op, n)$value, tolerance = 1e-11)
  }
})

test_that("several eigenvalues come back when the potential has several", {
  ## A deep well and a high barrier, well separated, give one value below the band
  ## and one above it.
  op <- finite_section(c(-6, 0, 0, 0, 6), at = -2:2)
  fit <- eigenpairs(op, n = 30L, k = 2L)
  expect_equal(nrow(fit$ritz), 2L)
  expect_true(all(fit$ritz$certified))
  expect_length(fit$values, 2L)
  expect_true(any(fit$values > 2))
  expect_true(any(fit$values < -2))
  ## each against the dense block it was computed from
  e <- eigen(explicit_block(op, 30L), symmetric = TRUE)$values
  expect_equal(sort(fit$values), sort(e[abs(e) > 2]), tolerance = 1e-10)
})

test_that("eigenpairs validates its arguments", {
  H <- finite_section(1)
  expect_error(eigenpairs(diag(3), n = 5), "expects a finite_section")
  expect_error(eigenpairs(H, n = 5, k = 0), "positive integer")
  expect_error(eigenpairs(H, n = 5, tol = 0), "positive number")
  expect_error(eigenpairs(H, n = 5, tol = c(1e-8, 1e-9)), "positive number")
  expect_error(eigenpairs(H, n = 0), "does not exceed the support radius")
})

test_that("print names what was certified and what was not", {
  expect_output(print(eigenpairs(finite_section(1), n = 12L)),
                "1 of 1 value certified")
  expect_output(print(eigenpairs(finite_section(1), n = 12L)),
                "dist to sigma\\(H\\) <=")
  expect_output(print(eigenpairs(finite_section(0), n = 12L)),
                "inside the band, not an eigenvalue")
})
