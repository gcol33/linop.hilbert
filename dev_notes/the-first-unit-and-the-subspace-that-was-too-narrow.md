# The first unit, and the subspace that was too narrow

Run 2026-08-02, against linop at `c7a82aa`. Everything below is measured; the probe
scripts are `dev_notes/spikes/`.

The unit is S0.6's class and it works: 319 assertions, `R CMD check` at one WARNING
whose two halves are both "linop is not on CRAN yet" and "the GitHub repo does not
exist yet". Three things came out of building it that were not in the design.

## 1. `linop::eigs()`'s default subspace is too narrow for this class

`eigs()` defaults to `ncv = max(2k + 1, k + 20)`, which is 21 at `k = 1`. That is a
reasonable default for a general hermitian operator and it is the wrong one here,
because an isolated eigenvalue of a finite section sits *just* outside a band of
essential spectrum and the block's remaining `2n` values cluster against it from
below. The eigenvalue is extremal but barely separated, which is the case a small
Krylov subspace handles worst.

`finite_section(0.3)`, whose eigenvalue is 2.022374841616 against a band edge of 2:

| n | ncv | error | eta | pairs converged |
|---|---|---|---|---|
| 40 | 21 | 4.241e-07 | 3.567e-04 | 1 of 1 |
| 80 | 21 | **7.347e-04** | 2.364e-03 | **0 of 1** |
| 80 | 40 | 2.724e-12 | 9.040e-07 | 1 of 1 |
| 80 | 80 | 2.724e-12 | 9.040e-07 | 1 of 1 |
| 120 | 21 | 4.441e-16 | 2.291e-09 | 1 of 1 |

**Widening the block made the answer three orders worse**, at `n = 80`, and the
only signal was `nconv`. So `eigenpairs()` defaults to `ncv = min(2n+1, max(40, 4k+20))`.
Forty was where it stopped mattering: `ncv = 80` and `ncv = 160` changed no digit
in any case tried, at either potential.

## 2. The certificate needed a row about the block, and it is a qualification

The failure above is invisible in the four rows the design called for. The bound is
still *true* -- that is the point -- so nothing can fail, and a caller reading a
number that got worse with more work would see a clean certificate.

The row is a qualification and never a failure, and the reason is worth keeping in
view. The residual bound

```
dist(theta, sigma(H)) <= ||(H - theta) u|| / ||u||
```

holds for **every** `theta` and every nonzero `u`, not only for eigenpairs. A
stalled block solve yields a true bound; it yields a worse one. What the row
reports is which of the two terms is binding, and `wider ncv is what moves it` is
the actionable half.

## 3. `eta` stops measuring truncation, and then raising `n` buys nothing

`eta / rho^n` is a constant of the potential while `eta` is measuring truncation.
For `finite_section(1)`, whose decay rate is the golden-ratio conjugate:

| n | eta | eta / rho^n |
|---|---|---|
| 20 | 3.864e-05 | 0.5845 |
| 40 | 2.554e-09 | 0.5844 |
| 80 | 7.734e-17 | 4.05 |
| 140 | 1.6e-13 | ~1e17 |

Past `n` around 60 the analytic residual has decayed below what the block eigensolve
resolves, and `eta` measures the eigensolve instead. It stops falling, and it stops
being monotone: 7.7e-17 at 80 and 1.6e-13 at 140, on the same operator with the same
settings. The value itself stays exact throughout, at 4.4e-16 from `sqrt(5)`.

So there is a width past which the reported bound cannot improve, and it is set by
the inner solve rather than by the mathematics. Documented on `eigenpairs()`'s `ncv`
argument, asserted in `test-eigenpairs.R`, and it is the reason the floor test
sweeps for the crossing instead of naming an `n`: **where the arithmetic floor bites
depends on how well the block was solved, so a test that fixed `n` would be testing
the eigensolver's convergence and not the floor.** With `ncv = 21` the crossing never
happens at all, because `eta` plateaus three orders above the true error.

At the default `ncv` it happens where S0.6 said: `n = 80`, `eta` 7.734e-17 against a
true error of 8.882e-16. The bound without the floor is violated on a fully
converged answer, which is the failure mode the floor exists for.

## What the unit needed from linop, and nothing else

`linop()`, `set_provenance()`, `provenance()`, the four provenance generics,
`eigs()`, `verify()`, `cap()`, `evidence()`, `requirement()`,
`evidence_satisfies()`, and the two names linop exported for exactly this:
`cert_rows()` and `build_certificate()`. No `:::` anywhere. The export decision was
scoped correctly.

## Corrections

| Where | Was | Is |
|---|---|---|
| `eigenpairs()` first draft | pass `ncv` through to `eigs()` and take its default | Its default converges nothing at `n = 80` on a weakly bound state. The default here is 40, measured |
| The certificate's four rows | floor, original residual, truncation bound, isolation gap | Five. A block eigensolve that stalled is invisible in those four and changes which term is binding |
| `test-eigenpairs.R` first draft | assert the floor bites at `n = 80` | It bites at a width that depends on `ncv`. The test sweeps and asserts the crossing happens somewhere, plus that the bound holds at every width |
