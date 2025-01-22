
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cheapr <a href="https://github.com/NicChr/cheapr"><img src="man/figures/cheapr_logo.png" align="right" height="238"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicChr/cheapr/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/cheapr)](https://CRAN.R-project.org/package=cheapr)
[![Codecov test
coverage](https://codecov.io/gh/NicChr/cheapr/graph/badge.svg)](https://app.codecov.io/gh/NicChr/cheapr)

<!-- badges: end -->

In cheapr, ‘cheap’ means fast and memory-efficient, and that’s exactly
the philosophy that cheapr aims to follow.

## Installation

You can install cheapr like so:

``` r
install.packages("cheapr")
```

or you can install the development version of cheapr:

``` r
remotes::install_github("NicChr/cheapr")
```

Some common operations that cheapr can do much faster and more
efficiently include:

- Counting, finding, removing and replacing `NA` and scalar values

- Creating factors

- Creating multiple sequences in a vectorised way

- Sub-setting vectors and data frames efficiently

- Safe, flexible and fast greatest common divisor and lowest common
  multiple

- Lags/leads

- Lightweight `integer64` support

- In-memory Math (no copies, vectors updated by reference)

- Summary statistics of data frame variables

- Binning of continuous data

Let’s first load the required packages

``` r
library(cheapr)
library(bench)
```

`num_na()` is a useful function to efficiently return the number of `NA`
values and can be used in a variety of problems.

Almost all the `NA` handling functions in cheapr can make use of
multiple cores on your machine through openMP.

``` r
x <- rep(NA, 10^6)

# 1 core by default
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)       1.69ms   1.75ms      555.    2.41KB      0  
#> 2 sum(is.na(x))  705.8µs   1.16ms      925.    3.81MB     77.3
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        464µs  475.5µs     2028.        0B      0  
#> 2 sum(is.na(x))    704µs   1.17ms      917.    3.81MB     75.5
options(cheapr.cores = 1)
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
#> 1 row_na_counts(m)    1.35ms   1.42ms      688.    9.16KB      0  
#> 2 rowSums(is.na(m))   2.44ms   2.91ms      350.    3.82MB     31.8
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)    1.53ms   1.58ms      604.    9.15KB      0  
#> 2 colSums(is.na(m))  996.6µs   1.47ms      705.    3.82MB     64.3
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6)
x[sample.int(10^6, 10^5)] <- NA
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      589µs   1.18ms      894.    3.81MB     84.1
#> 2 is_na(x)     1.02ms   1.72ms      614.    3.82MB     57.5

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
#> 1 is.na(hours)  997.8ms  997.8ms      1.00      61MB     1.00
#> 2 is_na(hours)   12.1ms   13.9ms     69.3     13.9MB     7.92
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
set.seed(42)
df <- data.frame(
  x = sample.int(100, 10^7, TRUE),
  y = factor_(sample(LETTERS, 10^7, TRUE)),
  z = rnorm(10^7)
)
overview(df)
#> obs: 10000000 
#> cols: 3 
#> 
#> ----- Numeric -----
#>   col  class n_missng p_complt n_unique      mean    p0   p25      p50   p75
#> 1   x integr        0        1      100     50.51     1    26       51    76
#> 2   z numerc        0        1 10000000 -0.000089 -5.16 -0.68 -0.00014  0.67
#>     p100   iqr    sd  hist
#> 1    100    50 28.86 ▇▇▇▇▇
#> 2   5.26  1.35     1 ▁▂▇▂▁
#> 
#> ----- Categorical -----
#>   col  class n_missng p_complt n_unique n_levels min max
#> 1   y factor        0        1       26       26   A   Z
mark(overview(df, hist = FALSE))
#> # A tibble: 1 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 overview(df, hist = FALSE)    1.32s    1.32s     0.757    2.09KB        0
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
#> 1 sset(x, x %in_% y)    108µs    160µs     6161.    83.2KB     2.22
#> 2 sset(x, x %in% y)     189µs    274µs     3443.   285.4KB     4.43
#> 3 x[x %in% y]           160µs    233µs     4278.   324.5KB     6.80
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            149.4µs  399.4µs     2658.    1.53MB    20.4 
#> 2 df[0:10^5, , drop = FALSE]    6.7ms   7.21ms      138.    4.83MB     2.06
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
#> 1 sset(df, -10^4:0)             26.9ms   34.3ms     24.8      152MB    15.3 
#> 2 df[-10^4:0, , drop = FALSE]  510.7ms  510.7ms      1.96     776MB     9.79
```

The biggest difference between `sset` and `[` is the way logical vectors
are handled. The two main differences when `i` is a logical vector are:

- `NA` values are ignored, only the locations of `TRUE` values are used.
- `i` must be the same length as `x` and is not recycled.

``` r
# Examples with NAs
x <- c(1, 5, NA, NA, -5)
x[x > 0]
#> [1]  1  5 NA NA
sset(x, x > 0)
#> [1] 1 5

# Example with length(i) < length(x)
sset(x, TRUE)
#> Error in check_length(i, length(x)): i must have length 5

# This is equivalent 
x[TRUE]
#> [1]  1  5 NA NA -5
# to..
sset(x)
#> [1]  1  5 NA NA -5
```

## Vector and data frame lags with `lag_()`

``` r
set.seed(37)
lag_(1:10, 3) # Lag(3)
#>  [1] NA NA NA  1  2  3  4  5  6  7
lag_(1:10, -3) # Lead(3)
#>  [1]  4  5  6  7  8  9 10 NA NA NA

# Using an example from data.table
library(data.table)
dt <- data.table(year=2010:2014, v1=runif(5), v2=1:5, v3=letters[1:5])

# Similar to data.table::shift()

lag_(dt, 1) # Lag 
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:    NA         NA    NA   <NA>
#> 2:  2010 0.54964085     1      a
#> 3:  2011 0.07883715     2      b
#> 4:  2012 0.64879698     3      c
#> 5:  2013 0.49685336     4      d
lag_(dt, -1) # Lead
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:  2011 0.07883715     2      b
#> 2:  2012 0.64879698     3      c
#> 3:  2013 0.49685336     4      d
#> 4:  2014 0.71878731     5      e
#> 5:    NA         NA    NA   <NA>
```

With `lag_` we can update variables by reference, including entire data
frames

``` r
# At the moment, shift() cannot do this
lag_(dt, set = TRUE)
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:    NA         NA    NA   <NA>
#> 2:  2010 0.54964085     1      a
#> 3:  2011 0.07883715     2      b
#> 4:  2012 0.64879698     3      c
#> 5:  2013 0.49685336     4      d

dt # Was updated by reference
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:    NA         NA    NA   <NA>
#> 2:  2010 0.54964085     1      a
#> 3:  2011 0.07883715     2      b
#> 4:  2012 0.64879698     3      c
#> 5:  2013 0.49685336     4      d
```

`lag2_` is a more generalised variant that supports vectors of lags,
custom ordering and run lengths.

``` r
lag2_(dt, order = 5:1) # Reverse order lag (same as lead)
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:  2010 0.54964085     1      a
#> 2:  2011 0.07883715     2      b
#> 3:  2012 0.64879698     3      c
#> 4:  2013 0.49685336     4      d
#> 5:    NA         NA    NA   <NA>
lag2_(dt, -1) # Same as above
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:  2010 0.54964085     1      a
#> 2:  2011 0.07883715     2      b
#> 3:  2012 0.64879698     3      c
#> 4:  2013 0.49685336     4      d
#> 5:    NA         NA    NA   <NA>
lag2_(dt, c(1, -1)) # Alternating lead/lag
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:    NA         NA    NA   <NA>
#> 2:  2011 0.07883715     2      b
#> 3:  2010 0.54964085     1      a
#> 4:  2013 0.49685336     4      d
#> 5:  2012 0.64879698     3      c
lag2_(dt, c(-1, 0, 0, 0, 0)) # Lead e.g. only first row
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:  2010 0.54964085     1      a
#> 2:  2010 0.54964085     1      a
#> 3:  2011 0.07883715     2      b
#> 4:  2012 0.64879698     3      c
#> 5:  2013 0.49685336     4      d
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
#> 1 gcd(x)        800ns    900ns   955932.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       36.4ms   37.2ms      26.7        0B        0
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
sequence_(seq_sizes, from = 0, by = 1/3, add_id = TRUE)
#>         1         1         1         2         2         2         2         2 
#> 0.0000000 0.3333333 0.6666667 0.0000000 0.3333333 0.6666667 1.0000000 1.3333333 
#>         3         3         3         3         3         3         3         3 
#> 0.0000000 0.3333333 0.6666667 1.0000000 1.3333333 1.6666667 2.0000000 2.3333333 
#>         3         3 
#> 2.6666667 3.0000000
```

You can also calculate the sequence sizes using `seq_size()`

``` r
seq_size(start, end, increments)
#> [1]  5  9 41
```

## ‘Cheaper’ Base R alternatives

### which

``` r
x <- rep(TRUE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    5.3ms   5.64ms      165.    3.81MB      0  
#> 2 base_which    559.8µs    1.5ms      718.    7.63MB     14.1
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.41ms   1.74ms      572.        0B      0  
#> 2 base_which    224.6µs    242µs     3990.    3.81MB     33.8
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   4.27ms    4.5ms      219.    1.91MB      0  
#> 2 base_which    538.8µs   1.02ms      982.    7.63MB     18.6
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   4.91ms   5.23ms      188.    3.81MB      0  
#> 2 base_which    751.5µs   1.76ms      599.    9.54MB     15.8
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   3.28ms   3.58ms      275.    1.89MB     0   
#> 2 base_which     3.53ms   4.03ms      245.     5.7MB     4.46
```

### factor

``` r
x <- sample(seq(-10^3, 10^3, 0.01))
y <- do.call(paste0, expand.grid(letters, letters, letters, letters))
mark(cheapr_factor = factor_(x), 
     base_factor = factor(x))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   8.31ms   9.38ms    104.      4.59MB        0
#> 2 base_factor   284.17ms 289.45ms      3.45   27.84MB        0
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   3.28ms   3.97ms    252.      1.53MB        0
#> 2 base_factor   426.68ms 443.79ms      2.25   22.79MB        0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 189.76ms 198.89ms     5.09     5.23MB        0
#> 2 base_factor      2.58s    2.58s     0.387   54.35MB        0
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   5.16ms   6.37ms     156.     3.49MB     0   
#> 2 base_factor    48.76ms  51.87ms      19.2   39.89MB     2.75
```

### intersect & setdiff

``` r
x <- sample.int(10^6, 10^5, TRUE)
y <- sample.int(10^6, 10^5, TRUE)
mark(cheapr_intersect = intersect_(x, y, dups = FALSE),
     base_intersect = intersect(x, y))
#> # A tibble: 2 × 6
#>   expression            min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>       <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_intersect   2.13ms   2.21ms      439.    1.18MB     0   
#> 2 base_intersect     3.73ms   4.21ms      232.    5.16MB     2.28
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   2.27ms   2.35ms      402.    1.76MB     2.26
#> 2 base_setdiff     3.84ms   4.21ms      234.    5.71MB     2.25
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.38ms   1.45ms      642.  781.34KB     0   
#> 2 base         2.25ms    2.5ms      388.    2.53MB     2.25
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.31ms   1.43ms      689.  787.87KB     0   
#> 2 base         2.22ms   2.63ms      377.    2.91MB     2.23
```

### `as_discrete`

`as_discrete` is a cheaper alternative to `cut`

``` r
x <- rnorm(10^7)
b <- seq(0, max(x), 0.2)
mark(cheapr_cut = as_discrete(x, b, left = FALSE), 
     base_cut = cut(x, b))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_cut    201ms    202ms      4.96    38.2MB     0   
#> 2 base_cut      515ms    515ms      1.94   267.1MB     1.94
```

### `cheapr_if_else`

A cheap alternative to `ifelse`

``` r
mark(
  cheapr_if_else(x >= 0, "pos", "neg"),
  ifelse(x >= 0, "pos", "neg"),
  data.table::fifelse(x >= 0, "pos", "neg")
)
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 3 × 6
#>   expression                           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 "cheapr_if_else(x >= 0, \"pos\… 111.77ms 118.24ms     8.32      114MB    1.66 
#> 2 "ifelse(x >= 0, \"pos\", \"neg…    1.25s    1.25s     0.803     534MB    0.803
#> 3 "data.table::fifelse(x >= 0, \… 105.92ms 113.03ms     8.58      114MB    3.43
```

### `case`

cheapr’s version of a case-when statement, with mostly the same
arguments as `dplyr::case_when` but similar efficiency as
`data.table::fcase`

``` r
mark(case(
    x >= 0 ~ "pos", 
    x < 0 ~ "neg", 
    .default = "Unknown"
),
data.table::fcase(
    x >= 0, "pos", 
    x < 0, "neg", 
    rep_len(TRUE, length(x)), "Unknown"
))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 "case(x >= 0 ~ \"pos\", x < 0 ~ \"n… 353ms  367ms      2.72     286MB     2.72
#> 2 "data.table::fcase(x >= 0, \"pos\",… 191ms  221ms      4.37     267MB     2.91
```

`val_match` is an even cheaper special variant of `case` when all LHS
expressions are length-1 vectors, i.e scalars

``` r
x <- round(rnorm(10^7))

mark(
  val_match(x, 1 ~ Inf, 2 ~ -Inf, .default = NaN),
     case(x == 1 ~ Inf, 
          x == 2 ~ -Inf, 
          .default = NaN),
     data.table::fcase(x == 1, Inf, 
          x == 2, -Inf, 
          rep_len(TRUE, length(x)), NaN)
     )
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 3 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                        <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 val_match(x, 1 ~ Inf, 2 ~ -Inf, …  97.5ms  99.6ms      9.64    87.9MB     1.93
#> 2 case(x == 1 ~ Inf, x == 2 ~ -Inf… 312.7ms 316.3ms      3.16   276.2MB     3.16
#> 3 data.table::fcase(x == 1, Inf, x… 156.3ms 161.6ms      6.21   305.2MB     6.21
```

`get_breaks` is a very fast function for generating pretty equal-width
breaks It is similar to `base::pretty` though somewhat less flexible
with simpler arguments.

``` r
x <- with_local_seed(rnorm(10^5), 112)
# approximately 10 breaks
get_breaks(x, 10)
#> [1] -6 -4 -2  0  2  4  6
pretty(x, 10)
#>  [1] -6 -5 -4 -3 -2 -1  0  1  2  3  4  5

mark(
  get_breaks(x, 20),
  pretty(x, 20), 
  check = FALSE
)
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 get_breaks(x, 20)   63.2µs   68.4µs    13990.        0B     0   
#> 2 pretty(x, 20)      408.9µs    663µs     1534.    1.91MB     9.78

# Not pretty but equal width breaks
get_breaks(x, 5, pretty = FALSE)
#> [1] -5.0135893 -3.2004889 -1.3873886  0.4257118  2.2388121  4.0519125
diff(get_breaks(x, 5, pretty = FALSE)) # Widths
#> [1] 1.8131 1.8131 1.8131 1.8131 1.8131
```

It can accept both data and a length-two vector representing a range,
meaning it can easily be used in ggplot2 and base R plots

``` r
library(ggplot2)
gg <- airquality |> 
    ggplot(aes(x = Ozone, y = Wind)) +
    geom_point() + 
    geom_smooth(se = FALSE)

# Add our breaks
gg +
  scale_x_continuous(breaks = get_breaks)
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
#> Warning: Removed 37 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
#> Warning: Removed 37 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```

<img src="man/figures/README-unnamed-chunk-30-1.png" width="100%" />

``` r

# More breaks

# get_breaks accepts a range too
gg +
  scale_x_continuous(breaks = \(x) get_breaks(range(x), 20)) 
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
#> Warning: Removed 37 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
#> Removed 37 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```

<img src="man/figures/README-unnamed-chunk-30-2.png" width="100%" />
