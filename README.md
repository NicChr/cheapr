
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cheapr

<!-- badges: start -->

[![R-CMD-check](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/cheapr)](https://CRAN.R-project.org/package=cheapr)
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
  # num_na is recursive so we compare it to unlisted length
  if (num_na(x) %in% c(0, unlisted_length(x))){
    x
  } else {
    vec_fill_missing(x, direction = "down")
  }
}
x <- rep(NA, 10^6)
identical(x, na_locf(x))
#> [1] TRUE
mark(na_locf(x), vec_fill_missing(x, direction = "down"))
#> # A tibble: 2 × 6
#>   expression                           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 "na_locf(x)"                     841.5µs  853.7µs     1156.        0B       0 
#> 2 "vec_fill_missing(x, direction…   2.67ms   2.84ms      346.    11.4MB     118.
mark(na_locf(x), vec_fill_missing(x, direction = "down"))
#> # A tibble: 2 × 6
#>   expression                           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 "na_locf(x)"                     841.2µs 847.15µs     1162.        0B       0 
#> 2 "vec_fill_missing(x, direction…   2.61ms   2.76ms      356.    11.4MB     185.
```

All the `NA` handling functions in cheapr can make use of multiple cores
on your machine using openMP.

``` r
# 1 core by default
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        839µs  844.7µs     1175.        0B      0  
#> 2 sum(is.na(x))    911µs   1.07ms      921.    3.81MB     84.7
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        242µs  305.4µs     2950.        0B      0  
#> 2 sum(is.na(x))    908µs   1.08ms      905.    3.81MB     82.3
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
#> 1 row_na_counts(m)    1.61ms    3.3ms      309.    12.9KB      0  
#> 2 rowSums(is.na(m))   2.79ms   2.88ms      344.    3.82MB     35.2
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)   693.7µs  846.4µs     1195.    12.9KB      0  
#> 2 colSums(is.na(m))   1.92ms   2.08ms      472.    3.82MB     48.1
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6)
x[sample.int(10^6, 10^5)] <- NA
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)     1.04ms   1.09ms      900.    3.81MB     145.
#> 2 is_na(x)    584.9µs  732.8µs     1345.    3.82MB     193.

### posixlt method is much faster
hours <- as.POSIXlt(seq.int(0, length.out = 10^6, by = 3600),
                    tz = "UTC")
hours[sample.int(10^6, 10^5)] <- NA

mark(is.na(hours), is_na(hours))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(hours)    1.22s    1.22s     0.821   61.05MB     1.64
#> 2 is_na(hours)   5.43ms   5.81ms   162.       3.83MB     7.99
```

It differs in 2 regards:

- List elements are regarded as `NA` when either that element is an `NA`
  value or it is a list containing only `NA` values.
- For data frames, `is_na` returns a logical vector where `TRUE` defines
  an empty row of only `NA` values.

``` r
# List example
is.na(list(NA, list(NA, NA), 10))
#> [1]  TRUE FALSE FALSE
is_na(list(NA, list(NA, NA), 10))
#> [1]  TRUE  TRUE FALSE

# Data frame example
df <- data.frame(x = c(1, NA, 3),
                 y = c(NA, NA, NA))
df
#>    x  y
#> 1  1 NA
#> 2 NA NA
#> 3  3 NA
is_na(df)
#> [1] FALSE  TRUE FALSE
is_na(df)
#> [1] FALSE  TRUE FALSE
# The below identity should hold
identical(is_na(df), row_na_counts(df) == ncol(df))
#> [1] TRUE
```

`is_na` and all the `NA` handling functions fall back on calling
`is.na()` if no suitable method is found. This means that custom objects
like vctrs rcrds and more are supported.

## Cheap data frame summaries with `overview`

Inspired by the excellent skimr package, `overview()` is a cheaper
alternative designed for larger data.

``` r
df <- data.frame(
  x = sample.int(100, 10^7, TRUE),
  y = factor_(sample(LETTERS, 10^7, TRUE)),
  z = rnorm(10^7)
)
overview(df, hist = TRUE)
#> rows: 10000000 cols: 3 
#> 
#> ----- Numeric -----
#>   col   class n_missing p_complete n_unique  mean    p0   p25 p50  p75 p100
#> 1   x integer         0          1      100 50.51     1    26  51   76  100
#> 3   z numeric         0          1 10000000     0 -5.29 -0.67   0 0.67 5.39
#>    iqr    sd  hist
#> 1   50 28.87 ▇▇▇▇▇
#> 3 1.35     1 ▁▂▇▂▁
#> 
#> ----- Categorical -----
#>   col  class n_missing p_complete n_unique n_levels min max
#> 2   y factor         0          1       26       26   A   Z
mark(overview(df))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 1 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 overview(df)    950ms    950ms      1.05    76.3MB     1.05
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
#> 1 gcd(x)        1.1µs    1.3µs   695642.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       54.7ms   54.8ms      18.3        0B        0
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
#> 1 cheapr_which   2.78ms   3.06ms      319.    3.82MB     6.47
#> 2 base_which     1.14ms   1.21ms      806.    7.63MB    34.8
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    207µs    271µs     3382.        0B      0  
#> 2 base_which      456µs    466µs     2120.    3.81MB     42.9
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.61ms   1.75ms      552.    1.91MB     6.35
#> 2 base_which     1.02ms   1.06ms      922.    7.63MB    34.8
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    3.7ms   3.82ms      256.    3.81MB     6.45
#> 2 base_which     1.36ms   1.44ms      671.    9.54MB    35.0
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   2.33ms   2.44ms      400.    1.89MB     2.05
#> 2 base_which     3.33ms   3.36ms      295.     5.7MB    11.1
```

### factor

``` r
# factor()
x <- sample(seq(-10^3, 10^3, 0.01))
y <- do.call(paste0, expand.grid(letters, letters, letters, letters))
mark(cheapr_factor = factor_(x), 
     base_factor = factor(x))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   10.1ms   10.6ms     94.1     4.59MB     2.14
#> 2 base_factor    598.5ms  598.5ms      1.67   27.84MB     0
mark(base_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 base_factor   5.56ms   5.87ms    171.      1.53MB     2.16
#> 2 base_factor 922.19ms 922.19ms      1.08   22.79MB     0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 195.03ms 195.31ms     5.09     5.23MB    0    
#> 2 base_factor      2.81s    2.81s     0.356   54.35MB    0.356
mark(base_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 base_factor    7.5ms   7.96ms     125.     3.49MB     2.11
#> 2 base_factor   43.1ms  45.55ms      22.0   39.89MB     9.43
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
#> 1 cheapr_intersect   3.03ms   3.36ms      296.    1.19MB     2.11
#> 2 base_intersect     4.38ms   4.61ms      215.    5.16MB     7.08
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   3.38ms    3.6ms      271.    1.76MB     4.48
#> 2 base_setdiff      4.7ms   4.93ms      202.    5.71MB     7.22
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.87ms   2.02ms      485.  785.44KB     4.47
#> 2 base         2.59ms   2.81ms      350.    2.53MB     4.49
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.87ms   2.01ms      491.  787.85KB     4.51
#> 2 base         2.68ms   2.99ms      333.    2.91MB     4.50
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
#> 1 cheapr_cut    131ms    131ms      7.61    38.1MB     2.54
#> 2 base_cut      409ms    409ms      2.45   267.1MB     2.45
```
