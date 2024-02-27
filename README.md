
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cheapr

<!-- badges: start -->
<!-- badges: end -->

In cheapr, ‘cheap’ means fast and memory-efficient, and that’s exactly
the philosophy that cheapr aims to follow.

## Installation

You can install the development version of cheapr like so:

``` r
remotes::install_github("NicChr/cheapr")
```

## Last-observation carried forward (minor optimisation)

Here is an example of a minor optimisation we can add to
`vctrs::vec_fill_missing` to return x if x has zero or only `NA` values.

``` r
library(cheapr)
library(vctrs)
#> Warning: package 'vctrs' was built under R version 4.3.2
library(bench)

na_locf <- function(x){
  if (num_na(x) %in% c(0, vec_size(x))){
    x
  } else {
    vec_fill_missing(x, direction = "down")
  }
}

x <- rep(NA, 10^6)
mark(na_locf(x), vec_fill_missing(x, direction = "down"))
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 "na_locf(x)"                        2.02ms 2.04ms      441.    4.68KB      0  
#> 2 "vec_fill_missing(x, direction = \…  5.6ms 5.99ms      163.   11.45MB     55.5
```

All the `NA` handling functions in cheapr can make use of multiple cores
on your machine using openMP.

``` r
# 1 core by default
mark(num_na(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)    2.01ms   2.01ms      458.        0B        0
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)     961µs   1.11ms      879.        0B        0
```

## Efficient NA counts by row/col

``` r
m <- matrix(x, ncol = 10^3)
# Number of NA values by row
mark(row_na_counts(m), 
     rowSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 row_na_counts(m)    3.48ms   6.74ms      153.    12.9KB      0  
#> 2 rowSums(is.na(m))   6.55ms   6.84ms      139.    3.82MB     15.8
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)    1.63ms   2.84ms      404.   21.02KB      0  
#> 2 colSums(is.na(m))   4.31ms   4.64ms      213.    3.82MB     22.1
```

## Greatest common divisor and smallest common multiple

``` r
gcd2(5, 25)
#> [1] 5
scm2(5, 6)
#> [1] 30

gcd(seq(5, 25, by = 5))
#> [1] 5
scm(seq(5, 25, by = 5))
#> [1] 300

x <- seq(1L, 1000000L, 1L)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)          3µs    3.4µs   271383.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)        112ms    116ms      8.43        0B        0
```

## ‘Cheaper’ Base R alternatives

### which

``` r
# which()
x <- rep(TRUE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   6.33ms   6.88ms      137.    3.82MB     21.5
#> 2 base_which     2.25ms   2.33ms      410.    7.63MB    135.
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which  502.8µs 904.75µs     1183.        0B      0  
#> 2 base_which     1.09ms   1.09ms      847.    3.81MB     79.3
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   4.09ms   4.72ms      212.    1.91MB     11.9
#> 2 base_which     2.23ms   2.34ms      404.    7.63MB    102.
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   8.77ms   9.29ms      104.    3.81MB     28.5
#> 2 base_which     2.79ms   2.88ms      327.    9.54MB    113.
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   5.75ms   6.12ms      163.    1.89MB     6.42
#> 2 base_which     7.84ms   8.07ms      121.     5.7MB    18.4
```

### factor

``` r
# factor()
x <- sample(seq(-10^3, 10^3, 0.01))
y <- do.call(paste0, expand.grid(letters, letters, letters, letters))
mark(cheapr_factor = factor_(x), 
     base_factor = factor(x))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor  20.32ms  22.89ms    38.9      10.5MB    3.89 
#> 2 base_factor      1.41s    1.41s     0.709    27.8MB    0.709
mark(base_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 base_factor  11.05ms  11.74ms    83.1      1.53MB     5.04
#> 2 base_factor    2.15s    2.15s     0.465   22.79MB     0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 394.39ms 415.26ms     2.41     9.23MB    1.20 
#> 2 base_factor      6.69s    6.69s     0.149   54.35MB    0.448
mark(base_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 base_factor   14.7ms     16ms     49.9     3.49MB     5.98
#> 2 base_factor  122.4ms    142ms      5.78   39.89MB    11.6
```

### intersect & setdiff

``` r
# intersect() & setdiff()
x <- sample.int(10^6, 10^5, TRUE)
y <- sample.int(10^6, 10^5, TRUE)
mark(cheapr_intersect = intersect_(x, y, dups = FALSE),
     base_intersect = intersect(x, y))
#> # A tibble: 2 × 6
#>   expression            min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>       <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_intersect   6.72ms   7.17ms     135.     1.19MB     2.10
#> 2 base_intersect    10.09ms  10.57ms      94.1    5.16MB    10.8
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff    7.2ms   7.71ms     118.     1.78MB     4.61
#> 2 base_setdiff     10.6ms   11.4ms      84.5    5.71MB    13.6
```

### cut.default

``` r
# cut.default()
x <- rnorm(10^7)
b <- seq(0, max(x), 0.2)
mark(cheapr_cut = cut_numeric(x, b), 
     base_cut = cut(x, b))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_cut    313ms    313ms      3.19    38.2MB     3.19
#> 2 base_cut      933ms    933ms      1.07   267.1MB     0
```
