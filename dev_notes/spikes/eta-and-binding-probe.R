setwd("C:/Users/Gilles Colling/Documents/dev/linop.hilbert")
devtools::load_all(".", quiet = TRUE)
options(digits = 12)

cat("=== where is the floor load-bearing? V = delta_0, truth sqrt(5) ===\n")
truth <- sqrt(5)
cat(sprintf("%5s %14s %12s %12s %12s %6s\n", "n", "theta", "err", "eta", "floor", "eta<err"))
for (n in c(40, 60, 80, 100, 120, 160, 200)) {
  fit <- eigenpairs(finite_section(1), n = n)
  err <- abs(fit$values[1] - truth)
  eta <- fit$ritz$eta[1]
  fl <- fit$certificate$values$floor
  cat(sprintf("%5d %14.12f %12.3e %12.3e %12.3e %6s\n", n, fit$values[1], err, eta, fl, eta < err))
}

cat("\n=== weak potential: where does it bind? ===\n")
for (v in c(0.3, 0.5)) {
  cat(sprintf(" v = %g, lambda = %.9f, rho = %.6f\n", v, sqrt(v^2 + 4),
              (-v + sqrt(v^2 + 4)) / 2))
  for (n in c(5, 8, 12, 20, 40, 80, 150)) {
    fit <- eigenpairs(finite_section(v), n = n)
    cat(sprintf("   n=%4d  theta=%.9f  gap=%+.3e  certified=%s  err=%.3e\n",
                n, fit$ritz$value[1], fit$ritz$gap[1], fit$ritz$certified[1],
                abs(fit$ritz$value[1] - sqrt(v^2 + 4))))
  }
}
