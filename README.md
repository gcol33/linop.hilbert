# linop.hilbert

*certified finite sections of operators on sequence spaces*

[![R-CMD-check](https://github.com/gcol33/linop.hilbert/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/linop.hilbert/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Eigenvalues of an operator on an infinite sequence space, with a bound on how
far the answer is from the truth.**

An operator on `l^2(Z)` has no matrix. The usual move is to truncate it to a
finite block and diagonalise that, which gives a number and no idea how good the
number is. This package computes that block, computes its eigenvalues, and returns
each one with a certificate bounding its distance to the spectrum of the operator
you started from.

```r
H <- finite_section(1)          # (Hu)_j = u_{j-1} + u_{j+1} + V_j u_j, V = delta_0
eigenpairs(H, n = 20)
#> <finite_section_eigen> n = 20, block of 41 sites
#>   1 of 1 value certified outside the essential spectrum [-2, 2]
#>     +2.23606797601   dist to sigma(H) <= 3.864e-05   in [2.23606797601, 2.23606798233]
#>   certificate: qualified
```

The value is `sqrt(5) = 2.2360679775`, it sits inside the bracket, and the bound
holds. Widen the block and both tighten:

```r
eigenpairs(H, n = 40)$certificate
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         construction     identity              1
#> block eigensolve               pass         computation      identity              1
#> original residual              pass         computation      identity              1
#> truncation bound               pass         theorem          deterministic_bound   1
#> isolation gap                  pass         computation      identity              1
#> ------------------------------------------------------------------------------------------
#> overall                        pass
```

## What the bound is made of

For this class of operator all three ingredients are closed form.

**The essential spectrum is known before anything runs.** A finitely supported
potential is a compact perturbation of the free operator, so by Weyl's theorem
`sigma_ess(H) = [-2, 2]` whatever `V` is. Everything outside that interval is an
isolated eigenvalue.

**The truncation residual is two numbers.** Extend the block eigenvector by zero
and apply `H`. Every site inside the block cancels, every site past `n+1` is zero
term by term, and what survives sits on exactly two sites: `u_n` and `u_{-n}`. So
the residual against the original operator costs no apply, no norm estimate, and
nothing that grows with `n`. Self-adjointness turns it into a spectral bound.

**The bracket is quadratic.** With the band edge known rather than estimated,
Kato-Temple gives `eta^2 / (q - a)` where the linear bound gives `eta` — four
orders of magnitude at `n = 20`.

## The case where it refuses

`V = 0` is the free operator. Its spectrum is `[-2, 2]` and purely continuous: it
has no eigenvalues at all. Every finite block still has `2n+1` of them.

```r
eigenpairs(finite_section(0), n = 20)
#> <finite_section_eigen> n = 20, block of 41 sites
#>   0 of 1 value certified outside the essential spectrum [-2, 2]
#>     +1.99440759436   inside the band, not an eigenvalue of H
#>   certificate: fail
```

The value is returned, in `$ritz`, because it is a legitimate approximation to a
point of continuous spectrum and dropping it would discard a decision that belongs
to you. It stays out of `$values`, where it would read as an eigenvalue of `H`.
The same inequality that decides this, `|theta| - 2 > 0`, is the one Kato-Temple
needs in its denominator.

## Installation

```r
# install.packages("pak")
pak::pak("gcol33/linop.hilbert")
```

Requires [linop](https://github.com/gcol33/linop), which supplies the matrix-free
operator, the eigensolver and the certificate object.

## Scope

One operator class: self-adjoint Jacobi with a finitely supported potential, the
class for which the truncation error is a closed form rather than an existence
statement. Wider spaces, closed operators with domains, forms, and Galerkin and
collocation discretisations are not here.
