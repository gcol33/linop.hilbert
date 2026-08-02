test_that("the constructor centres the support and records what is known a priori", {
  H <- finite_section(1)
  expect_s3_class(H, "finite_section")
  expect_equal(H$at, 0L)
  expect_equal(H$support_radius, 0L)
  ## sigma_ess(H) = [-2, 2] for every V, by Weyl against the free operator
  expect_equal(H$band_edge, 2)
  ## ||H|| <= ||free|| + ||V||, and both terms are exact for this class
  expect_equal(H$norm_bound, 3)

  H3 <- finite_section(c(1.5, -0.5, 2))
  expect_equal(H3$at, -1:1)
  expect_equal(H3$support_radius, 1L)
  expect_equal(H3$norm_bound, 4)

  ## an even number of sites has no centre; the convention is documented by the
  ## default and asserted here so it cannot drift
  expect_equal(finite_section(c(1, 2))$at, 0:1)
})

test_that("sites are sorted with their potential, however they arrive", {
  H <- finite_section(c(2, 1.5, -0.5), at = c(1, -1, 0))
  expect_equal(H$at, -1:1)
  expect_equal(H$V, c(1.5, -0.5, 2))
})

test_that("the constructor refuses what it cannot certify", {
  expect_error(finite_section(numeric(0)), "non-empty")
  expect_error(finite_section(c(1, NA)), "finite")
  expect_error(finite_section(c(1, Inf)), "finite")
  expect_error(finite_section("x"), "numeric")
  expect_error(finite_section(1, at = c(0, 1)), "one per entry")
  expect_error(finite_section(c(1, 2), at = c(0.5, 1)), "integer sites")
  expect_error(finite_section(c(1, 2), at = c(0, 0)), "repeats the site")
})

test_that("a complex potential is refused by name, not coerced", {
  ## H is self-adjoint because the diagonal is real. A complex V would break that
  ## silently: every bound below rests on self-adjointness and nothing downstream
  ## would notice.
  expect_error(finite_section(complex(real = 1, imaginary = 1)), "self-adjoint")
})

test_that("V = 0 constructs, because the free operator is a legitimate operator", {
  ## It has no eigenvalues, which is a fact eigenpairs() reports rather than a
  ## reason to refuse the object.
  H <- finite_section(0)
  expect_s3_class(H, "finite_section")
  expect_equal(H$norm_bound, 2)
})

test_that("print says what the operator is without running anything", {
  expect_output(print(finite_section(1)), "l\\^2\\(Z\\)")
  expect_output(print(finite_section(1)), "essential spectrum")
  expect_output(print(finite_section(c(1, 2, 3))), "3 sites")
  expect_output(print(finite_section(1)), "1 site,")
})
