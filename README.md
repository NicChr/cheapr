
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cheapr

<!-- badges: start -->

[![R-CMD-check](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

In cheapr, ‘cheap’ means fast and memory-efficient, and that’s exactly
the philosophy that cheapr aims to follow.

## Installation

You can install the development version of cheapr like so:

``` r
remotes::install_github("NicChr/cheapr")
```

## Last-observation carried forward (minor optimisation)

`num_na()` is a useful function to efficiently return the number of `NA`
values and can be used in a variety of problems.

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
#>   expression                           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 "na_locf(x)"                     843.7µs 849.65µs     1155.    4.68KB       0 
#> 2 "vec_fill_missing(x, direction…   2.65ms   2.76ms      355.   11.45MB     121.
```

All the `NA` handling functions in cheapr can make use of multiple cores
on your machine using openMP.

``` r
# 1 core by default
mark(num_na(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)     838µs    839µs     1167.        0B        0
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)     232µs    293µs     3161.        0B        0
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
#> 1 row_na_counts(m)    2.26ms   3.16ms      317.    12.9KB      0  
#> 2 rowSums(is.na(m))   2.76ms   2.85ms      343.    3.82MB     34.3
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)     683µs 723.25µs     1295.   21.02KB      0  
#> 2 colSums(is.na(m))   1.94ms   2.04ms      488.    3.82MB     49.3
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
#> 1 gcd(x)        1.2µs    1.3µs   717381.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       46.1ms   46.4ms      21.4        0B        0
```

## Creating many sequences

As an example, to create 3 sequences with different increments,  
the usual approach might be to use lapply to loop through the increment
values together with `seq()`

``` r
# Base R
increments <- c(1, 0.5, 0.1)
start <- 1
end <- 5
unlist(lapply(increments, \(x) seq(start, end, x)))
#>  [1] 1.0 2.0 3.0 4.0 5.0 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 1.0 1.1 1.2 1.3 1.4
#> [20] 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3
#> [39] 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0
```

In cheapr you can use `seq_()` which accepts vector arguments.

``` r
seq_(start, end, increments)
#>  [1] 1.0 2.0 3.0 4.0 5.0 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 1.0 1.1 1.2 1.3 1.4
#> [20] 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3
#> [39] 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0
```

Use `add_id = TRUE` to label the individual sequences.

``` r
seq_(start, end, increments, add_id = TRUE)
#>   1   1   1   1   1   2   2   2   2   2   2   2   2   2   3   3   3   3   3   3 
#> 1.0 2.0 3.0 4.0 5.0 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 1.0 1.1 1.2 1.3 1.4 1.5 
#>   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3 
#> 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 
#>   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3 
#> 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0
```

If you know the sizes of your sequences beforehand, use `sequence_()`

``` r
seq_sizes <- c(3, 5, 10)
sequence_(seq_sizes, from = 0, by = 1/3, add_id = TRUE) |> 
  enframe_()
#> # A tibble: 18 × 2
#>    name  value
#>    <chr> <dbl>
#>  1 1     0    
#>  2 1     0.333
#>  3 1     0.667
#>  4 2     0    
#>  5 2     0.333
#>  6 2     0.667
#>  7 2     1    
#>  8 2     1.33 
#>  9 3     0    
#> 10 3     0.333
#> 11 3     0.667
#> 12 3     1    
#> 13 3     1.33 
#> 14 3     1.67 
#> 15 3     2    
#> 16 3     2.33 
#> 17 3     2.67 
#> 18 3     3
```

You can also calculate the sequence sizes using `seq_size()`

``` r
seq_size(start, end, increments)
#> [1]  5  9 41
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
#> 1 cheapr_which   2.86ms   3.08ms      312.    3.82MB     62.9
#> 2 base_which     1.09ms   1.15ms      853.    7.63MB    405.
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    202µs    257µs     3671.        0B       0 
#> 2 base_which      464µs    468µs     2092.    3.81MB     310.
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.61ms   1.79ms      536.    1.91MB     35.7
#> 2 base_which        1ms   1.05ms      933.    7.63MB    314.
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   3.71ms   3.94ms      251.    3.81MB     65.5
#> 2 base_which     1.35ms   1.43ms      686.    9.54MB    189.
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   2.38ms   2.46ms      399.    1.89MB     20.6
#> 2 base_which     3.31ms   3.35ms      296.     5.7MB     48.5
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
#> 1 cheapr_factor   9.58ms   10.2ms     79.4     4.63MB    11.9 
#> 2 base_factor   581.08ms  581.1ms      1.72   27.84MB     1.72
mark(base_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 base_factor   5.24ms   5.74ms    155.      1.53MB     5.95
#> 2 base_factor 924.62ms 924.62ms      1.08   22.79MB     1.08
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 193.66ms 197.72ms     5.08     9.23MB    0    
#> 2 base_factor      2.86s    2.86s     0.349   54.35MB    0.699
mark(base_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 base_factor    7.1ms   7.77ms     105.     3.49MB     9.86
#> 2 base_factor   61.4ms  65.12ms      11.8   39.89MB    13.8
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
#> 1 cheapr_intersect      3ms   3.31ms      296.    1.19MB     6.73
#> 2 base_intersect      4.2ms   4.58ms      219.    5.16MB    33.1
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   3.15ms    3.6ms      276.    1.78MB     9.45
#> 2 base_setdiff     4.12ms   4.89ms      199.    5.71MB    24.9
```

### `%in_` and `!%in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.89ms   1.99ms      488.  785.45KB     9.52
#> 2 base         2.36ms   2.76ms      352.    2.53MB    18.6
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr        1.8ms   1.95ms      507.  785.46KB     6.88
#> 2 base         2.38ms   2.94ms      335.    2.91MB    22.5
```

### cut.default

``` r
# cut.default()
x <- rnorm(10^7)
b <- seq(0, max(x), 0.2)
mark(cheapr_cut = cut_numeric(x, b), 
     base_cut = cut(x, b))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_cut    130ms    131ms      7.54    38.2MB     1.89
#> 2 base_cut      422ms    445ms      2.25   267.1MB     2.25
```
