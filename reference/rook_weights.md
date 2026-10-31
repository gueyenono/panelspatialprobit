# Rook contiguity weights on a regular grid

Convenience constructor used in examples, tests and simulations.

## Usage

``` r
rook_weights(nrow, ncol)
```

## Arguments

- nrow, ncol:

  Grid dimensions; units are numbered row by row.

## Value

A binary sparse `dgCMatrix` of dimension `nrow * ncol` with unit
identifiers `"1"`, `"2"`, ... as dimnames.

## Examples

``` r
rook_weights(2, 3)
#> 6 x 6 sparse Matrix of class "dgCMatrix"
#>   1 2 3 4 5 6
#> 1 . 1 . 1 . .
#> 2 1 . 1 . 1 .
#> 3 . 1 . . . 1
#> 4 1 . . . 1 .
#> 5 . 1 . 1 . 1
#> 6 . . 1 . 1 .
```
