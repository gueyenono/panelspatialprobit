# Prepare a spatial weights matrix

Coerces `W` to a sparse `dgCMatrix`, validates it and optionally
row-standardises it. Unit identifiers are taken from `ids`, or from the
row names of `W`, and are used by
[`psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
to match rows of `W` to the units in the data.

## Usage

``` r
psp_weights(W, ids = NULL, style = c("asis", "W"))
```

## Arguments

- W:

  A square matrix, a `Matrix` sparse matrix, an `spdep` `listw` or `nb`
  object.

- ids:

  Optional vector of unit identifiers, one per row of `W`.

- style:

  `"asis"` (default) keeps the weights unchanged; `"W"` row-standardises
  them (rows of isolates stay zero).

## Value

A `dgCMatrix` with unit identifiers as dimnames.

## Examples

``` r
W <- psp_weights(rook_weights(3, 3), style = "W")
Matrix::rowSums(W)
#> 1 2 3 4 5 6 7 8 9 
#> 1 1 1 1 1 1 1 1 1 
```
