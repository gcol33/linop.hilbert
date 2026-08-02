setwd("C:/Users/Gilles Colling/Documents/dev/linop.hilbert")
devtools::load_all(".", quiet = TRUE)
options(digits = 12)

cat("=== does ncv explain the eta plateau and the n=80 regression? ===\n")
for (v in c(1, 0.3)) {
  truth <- sqrt(v^2 + 4)
  cat(sprintf("\n v = %g, truth = %.12f\n", v, truth))
  cat(sprintf("%5s %5s %12s %12s %8s %6s\n", "n", "ncv", "err", "eta", "iters", "nconv"))
  for (n in c(40, 80, 120)) {
    for (ncv in c(21, 40, 80, 160)) {
      if (ncv > 2 * n + 1) next
      fit <- eigenpairs(finite_section(v), n = n, ncv = ncv)
      cat(sprintf("%5d %5d %12.3e %12.3e %8d %6d\n", n, ncv,
                  abs(fit$ritz$value[1] - truth), fit$ritz$eta[1],
                  fit$inner$iterations, fit$inner$nconv))
    }
  }
}

cat("\n=== with a wide subspace, does eta fall below err? ===\n")
truth <- sqrt(5)
for (n in c(60, 80, 100, 140)) {
  fit <- eigenpairs(finite_section(1), n = n, ncv = min(2 * n + 1, 120))
  err <- abs(fit$values[1] - truth)
  eta <- fit$ritz$eta[1]
  fl <- fit$certificate$values$floor
  cat(sprintf(" n=%4d err=%.3e eta=%.3e floor=%.3e  eta<err: %s  err<=eta: %s\n",
              n, err, eta, fl, eta < err, err <= eta))
}
