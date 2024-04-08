
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cheapr

<!-- badges: start -->

[![R-CMD-check](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/cheapr)](https://CRAN.R-project.org/package=cheapr)
[![Codecov test
coverage](https://codecov.io/gh/NicChr/cheapr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/NicChr/cheapr?branch=main)
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
#> 1 "na_locf(x)"                     841.5µs 851.15µs     1152.        0B       0 
#> 2 "vec_fill_missing(x, direction…   2.59ms   2.79ms      353.    11.4MB     117.
mark(na_locf(x), vec_fill_missing(x, direction = "down"))
#> # A tibble: 2 × 6
#>   expression                           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 "na_locf(x)"                     841.7µs  851.5µs     1164.        0B       0 
#> 2 "vec_fill_missing(x, direction…   2.62ms   2.77ms      353.    11.4MB     174.
```

All the `NA` handling functions in cheapr can make use of multiple cores
on your machine using openMP.

``` r
# 1 core by default
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        839µs  844.6µs     1178.        0B      0  
#> 2 sum(is.na(x))    975µs   1.08ms      909.    3.81MB     83.1
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        234µs  297.5µs     3099.        0B      0  
#> 2 sum(is.na(x))    992µs   1.08ms      885.    3.81MB     80.2
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
#> 1 row_na_counts(m)    1.52ms   3.29ms      319.    12.9KB      0  
#> 2 rowSums(is.na(m))    2.8ms   2.89ms      338.    3.82MB     34.6
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)   684.1µs  799.8µs     1222.    12.9KB      0  
#> 2 colSums(is.na(m))   1.99ms   2.06ms      479.    3.82MB     47.4
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6)
x[sample.int(10^6, 10^5)] <- NA
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)     1.05ms   1.09ms      902.    3.81MB     83.6
#> 2 is_na(x)    573.7µs  702.2µs     1323.    3.82MB    121.

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
#> 1 is.na(hours)    1.17s    1.17s     0.854      61MB    0.854
#> 2 is_na(hours)   5.24ms    5.7ms   163.        9.8MB    9.93
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
#> obs: 10000000 
#> cols: 3 
#> 
#> ----- Numeric -----
#>   col   class n_missing p_complete n_unique  mean    p0   p25 p50  p75 p100
#> 1   x integer         0          1      100 50.49     1    25  50   75  100
#> 2   z numeric         0          1 10000000     0 -5.03 -0.67   0 0.67 5.35
#>    iqr    sd  hist
#> 1   50 28.87 ▇▇▇▇▇
#> 2 1.35     1 ▁▂▇▁▁
#> 
#> ----- Categorical -----
#>   col  class n_missing p_complete n_unique n_levels min max
#> 1   y factor         0          1       26       26   A   Z
mark(overview(df))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 1 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 overview(df)    922ms    922ms      1.08    76.3MB     1.08
```

## Cheaper and consistent subsetting with `sset`

``` r
sset(iris, 1:5)
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 1          5.1         3.5          1.4         0.2  setosa
#> 2          4.9         3.0          1.4         0.2  setosa
#> 3          4.7         3.2          1.3         0.2  setosa
#> 4          4.6         3.1          1.5         0.2  setosa
#> 5          5.0         3.6          1.4         0.2  setosa
sset(iris, 1:5, j = "Species")
#>   Species
#> 1  setosa
#> 2  setosa
#> 3  setosa
#> 4  setosa
#> 5  setosa

# sset always returns a data frame when input is a data frame

sset(iris, 1, 1) # data frame
#>   Sepal.Length
#> 1          5.1
iris[1, 1] # not a data frame
#> [1] 5.1

x <- sample.int(10^6, 10^4, TRUE)
y <- sample.int(10^6, 10^4, TRUE)
mark(sset(x, x %in_% y), sset(x, x %in% y), x[x %in% y])
#> # A tibble: 3 × 6
#>   expression              min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>         <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(x, x %in_% y)   90.6µs    114µs     8443.    83.1KB     2.06
#> 2 sset(x, x %in% y)   151.2µs    220µs     4454.   285.3KB     6.74
#> 3 x[x %in% y]         158.4µs    205µs     4826.   324.4KB     6.78
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            365.3µs  451.8µs     2148.    1.53MB    15.2 
#> 2 df[0:10^5, , drop = FALSE]   6.83ms    7.1ms      141.    4.82MB     4.33
```

It also accepts negative indexes

``` r
mark(sset(df, -10^4:0), 
     df[-10^4:0, , drop = FALSE],
     check = FALSE) # The only difference is the row names
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression                       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, -10^4:0)             21.9ms   28.5ms     29.4      152MB    20.2 
#> 2 df[-10^4:0, , drop = FALSE]  653.8ms  653.8ms      1.53     776MB     7.65
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
#> 1 gcd(x)        1.2µs    1.3µs   667249.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       46.1ms   46.2ms      21.6        0B        0
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
#> 1 cheapr_which    2.5ms   2.72ms      353.    3.81MB     2.06
#> 2 base_which     1.15ms   1.24ms      747.    7.63MB    11.9
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    214µs    263µs     3667.        0B      0  
#> 2 base_which      457µs    461µs     2141.    3.81MB     17.7
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.45ms   1.63ms      595.    1.91MB     4.19
#> 2 base_which     1.02ms   1.08ms      911.    7.63MB    13.5
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   3.39ms   3.53ms      273.    3.81MB     4.30
#> 2 base_which     1.38ms   1.55ms      622.    9.54MB    11.5
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   2.09ms   2.21ms      438.    1.89MB     0   
#> 2 base_which     3.33ms   3.38ms      295.     5.7MB     4.24
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
#> 1 cheapr_factor   9.88ms   10.4ms     95.5     4.59MB        0
#> 2 base_factor   523.52ms  523.5ms      1.91   27.84MB        0
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor    5.6ms   5.92ms    167.      1.53MB        0
#> 2 base_factor    800.8ms 800.85ms      1.25   22.79MB        0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 202.93ms 203.91ms     4.89     5.23MB        0
#> 2 base_factor      2.84s    2.84s     0.352   54.35MB        0
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   7.06ms   8.07ms     121.     3.49MB     0   
#> 2 base_factor    44.53ms  48.78ms      20.5   39.89MB     2.28
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
#> 1 cheapr_intersect   3.16ms   3.43ms      290.    1.18MB     0   
#> 2 base_intersect     4.42ms   4.61ms      212.    5.16MB     2.21
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff    3.4ms   3.73ms      266.    1.76MB     0   
#> 2 base_setdiff     4.66ms      5ms      195.    5.71MB     2.22
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.96ms   2.05ms      481.  781.34KB     0   
#> 2 base         2.69ms   2.82ms      345.    2.53MB     2.25
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.95ms   2.09ms      476.  787.85KB     0   
#> 2 base         2.86ms   2.98ms      324.    2.91MB     2.23
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
#> 1 cheapr_cut    132ms    133ms      7.50    38.1MB     0   
#> 2 base_cut      438ms    486ms      2.06   267.1MB     2.06
```
