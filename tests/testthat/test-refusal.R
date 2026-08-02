## Section 8.9's case, and the one that decides whether this package is honest.
## The free operator has no eigenvalues at all: sigma(H_free) = [-2, 2] is purely
## continuous. Every finite section of it has 2n+1 eigenvalues, all real, all
## inside the band, converging to the edge from below. Nothing about the finite
## computation is wrong, and none of its output is an eigenvalue of H.

test_that("the free operator yields no certified value, at any width", {
  H <- finite_section(0)
  for (n in c(10L, 40L, 100L)) {
    fit <- eigenpairs(H, n = n)
    expect_length(fit$values, 0L)
    expect_false(fit$ritz$certified[1L])
    expect_lt(abs(fit$ritz$value[1L]), 2)
    expect_equal(fit$certificate$overall, "fail")
  }
})

test_that("the failing row is the isolation gap and it says why", {
  cert <- eigenpairs(finite_section(0), n = 20L)$certificate
  failed <- cert$checks[cert$checks$status == "fail", ]
  expect_equal(failed$check, "isolation gap")
  expect_match(failed$detail, "inside the essential spectrum")
  expect_match(failed$detail, "discretisations of\\s+continuous spectrum")
  ## the truncation bound is still computed and still true; what is missing is
  ## anything for it to be a bound about
  expect_true("truncation bound" %in% cert$checks$check)
})

test_that("the uncertified value is returned and is not called an eigenvalue", {
  ## Dropping it would destroy a decision the caller may want to make -- the value
  ## is a legitimate approximation to a point of continuous spectrum. Putting it in
  ## $values would assert something false. It goes in $ritz, labelled.
  fit <- eigenpairs(finite_section(0), n = 30L)
  expect_length(fit$values, 0L)
  expect_equal(nrow(fit$ritz), 1L)
  expect_false(fit$ritz$certified)
  expect_true(is.na(fit$ritz$lower))
  expect_true(is.na(fit$ritz$upper))
})

test_that("no error is thrown: a refusal is a certificate", {
  ## A contradicted declaration would be an error. This is not one -- the operator
  ## is exactly what it said it was, and the answer to the question asked is that
  ## there is no eigenvalue here.
  expect_no_error(eigenpairs(finite_section(0), n = 20L))
  expect_silent(invisible(eigenpairs(finite_section(0), n = 20L)))
})

test_that("a state that does not fit in the block is not certified at that width", {
  ## Every attractive well in one dimension binds, so the honest weak case is a
  ## state whose extent exceeds the block. A shallow well has a decay rate near 1,
  ## so at a narrow n the Ritz value is still inside the band and the answer is
  ## the same as for V = 0: not certified. The width is what changes it.
  v <- 0.3
  expect_gt(decay_rate(single_well_eigenvalue(v)), 0.85)  # slow decay, wide state

  narrow <- eigenpairs(finite_section(v), n = 5L)
  expect_false(narrow$ritz$certified[1L])
  expect_length(narrow$values, 0L)
  expect_lt(abs(narrow$ritz$value[1L]), 2)

  wide <- eigenpairs(finite_section(v), n = 120L)
  expect_true(wide$ritz$certified[1L])
  expect_lt(abs(wide$values[1L] - single_well_eigenvalue(v)), 1e-8)
})

test_that("asking for more values than the operator has certifies only the real ones", {
  ## One well, one bound state. The second value asked for is a band value, and the
  ## certificate fails on it while the first stays exact.
  fit <- eigenpairs(finite_section(3), n = 30L, k = 2L)
  expect_equal(sum(fit$ritz$certified), 1L)
  expect_length(fit$values, 1L)
  expect_lt(abs(fit$values[1L] - single_well_eigenvalue(3)), 1e-9)
  expect_equal(fit$certificate$overall, "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "isolation gap"],
               "1 of 2 values")
})
