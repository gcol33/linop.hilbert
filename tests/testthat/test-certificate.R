test_that("the certificate is the shape this package issues, not one of linop's", {
  cert <- eigenpairs(finite_section(1), n = 20L)$certificate
  expect_s3_class(cert, "linop_certificate")
  expect_equal(cert$subject, "finite section")
  expect_equal(cert$checks$check,
               c("arithmetic floor", "block eigensolve", "original residual",
                 "truncation bound", "isolation gap"))
  ## no residual row for the finite problem and no backward-error row: those are
  ## questions about H_n, and this certificate is about H
  expect_false("backward error" %in% cert$checks$check)
  expect_false("residual" %in% cert$checks$check)
})

test_that("the truncation bound carries the theorem and what the theorem rests on", {
  cert <- eigenpairs(finite_section(1), n = 20L)$certificate
  row <- cert$checks[cert$checks$check == "truncation bound", ]
  expect_equal(row$source, "theorem")
  expect_equal(row$guarantee, "deterministic_bound")

  ## the bound is a spectral statement only because H is self-adjoint, and for this
  ## class that is a property of the constructor rather than a declaration
  ev <- cert$evidence[["truncation bound"]]
  expect_s3_class(ev, "linop_evidence")
  expect_length(ev$depends_on, 1L)
  expect_equal(ev$depends_on[[1L]]$source, "construction")

  ## so it satisfies a requirement that refuses declarations, which the same bound
  ## resting on a user_declaration would not
  expect_true(linop::evidence_satisfies(
    ev, linop::requirement(sources = c("construction", "computation", "theorem"))))
})

test_that("the floor row states the norm exactly, since this class knows it", {
  op <- finite_section(c(1.5, -0.5, 2))
  cert <- eigenpairs(op, n = 15L)$certificate
  row <- cert$checks[cert$checks$check == "arithmetic floor", ]
  expect_equal(row$source, "construction")
  expect_equal(row$guarantee, "identity")
  ## ||H|| <= 2 + max|V| = 4, with no estimate anywhere
  expect_match(row$detail, "2 \\+ max\\|V\\| = 4 exactly")
  expect_equal(cert$values$floor, 4 * op$norm_bound * .Machine$double.eps)
})

test_that("every row is deterministic, which is the point of this class", {
  cert <- eigenpairs(finite_section(1), n = 40L)$certificate
  expect_equal(cert$without_deterministic_bound, character(0))
  expect_equal(cert$overall, "pass")
})

test_that("a bound that has not reached tol is qualified rather than passed", {
  ## n = 5 leaves eta at about 5e-2, four orders above the default tol. Nothing is
  ## wrong; the answer is simply not yet good, and the certificate says which of the
  ## two it is.
  cert <- eigenpairs(finite_section(1), n = 5L)$certificate
  expect_equal(cert$overall, "qualified")
  expect_equal(cert$checks$status[cert$checks$check == "original residual"], "qualified")
  expect_equal(cert$checks$status[cert$checks$check == "isolation gap"], "pass")
})

test_that("a stalled block eigensolve qualifies the certificate and never fails it", {
  ## The residual bound holds for any theta and any nonzero u, so a block solve
  ## that stalled still gives a true bound. What it gives is a worse one, and the
  ## row says which of the two limits the number being read.
  fit <- eigenpairs(finite_section(0.3), n = 80L, ncv = 21L)
  cert <- fit$certificate
  row <- cert$checks[cert$checks$check == "block eigensolve", ]
  expect_equal(row$status, "qualified")
  expect_match(row$detail, "0 of 1 pairs")
  expect_match(row$detail, "wider ncv")
  expect_false("fail" %in% cert$checks$status)

  ## and the bound it reports is still true, which is the point
  expect_lte(abs(fit$ritz$value[1L] - single_well_eigenvalue(0.3)),
             fit$ritz$bound[1L])

  ## the default ncv is what makes the same width work
  ok <- eigenpairs(finite_section(0.3), n = 80L)
  expect_equal(ok$certificate$checks$status[
    ok$certificate$checks$check == "block eigensolve"], "pass")
  expect_lt(abs(ok$values[1L] - single_well_eigenvalue(0.3)), 1e-9)
})

test_that("the inner eigensolve's own certificate is carried, not restated", {
  fit <- eigenpairs(finite_section(1), n = 20L)
  inner <- fit$certificate$values$inner
  expect_s3_class(inner, "linop_certificate")
  expect_equal(inner$subject, "eigen")
  ## it answers a different question, about H_n, and is kept for a caller who wants
  ## to know how well the block itself was solved
  expect_true("backward error" %in% inner$checks$check)
  expect_identical(inner, fit$inner$certificate)
})

test_that("summary prints both the result and the certificate", {
  fit <- eigenpairs(finite_section(1), n = 12L)
  expect_output(summary(fit), "finite_section_eigen")
  expect_output(summary(fit), "truncation bound")
  expect_output(summary(fit), "overall")
})
