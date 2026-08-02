stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

is_scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x)

is_whole <- function(x) {
  is.numeric(x) && all(is.finite(x)) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
}

## The package name appears in error messages and in the provenance envelope, and
## a typo in either is silent. One constant.
PROVIDER <- "linop.hilbert"
