
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
#> 
#> Attaching package: 'cheapr'
#> The following objects are masked from 'package:base':
#> 
#>     round, trunc
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
#> 1 num_na(x)        983µs   1.01ms      973.    2.41KB      0  
#> 2 sum(is.na(x))    780µs   1.84ms      543.    3.81MB     44.8
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        258µs  346.5µs     2666.        0B      0  
#> 2 sum(is.na(x))    812µs   1.93ms      493.    3.81MB     40.4
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
#> 1 row_na_counts(m)    1.89ms   2.21ms      455.    9.14KB      0  
#> 2 rowSums(is.na(m))   2.63ms   3.71ms      267.    3.82MB     25.8
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)    2.49ms   2.55ms      384.    9.14KB      0  
#> 2 colSums(is.na(m))   1.77ms   2.85ms      353.    3.82MB     31.0
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6)
x[sample.int(10^6, 10^5)] <- NA
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      818µs   2.01ms      481.    3.81MB     61.6
#> 2 is_na(x)     1.33ms    2.6ms      391.    3.82MB     36.2

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
#> 1 is.na(hours)    1.19s    1.19s     0.843      61MB    0.843
#> 2 is_na(hours)  13.26ms  15.35ms    64.3      13.9MB    7.80
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
#>   col class n_mssn p_cmpl   n_uniq  mean    p0   p25 p50  p75 p100  iqr    sd
#> 1   x intgr      0      1      100 50.51     1    26  51   76  100   50 28.86
#> 2   z numrc      0      1 10000000     0 -5.16 -0.68   0 0.67 5.26 1.35     1
#>    hist
#> 1 ▇▇▇▇▇
#> 2 ▁▂▇▂▁
#> 
#> ----- Categorical -----
#>   col class n_mssn p_cmpl n_uniq n_lvls min max
#> 1   y factr      0      1     26     26   A   Z
mark(overview(df, hist = FALSE))
#> # A tibble: 1 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 overview(df, hist = FALSE)    1.33s    1.33s     0.754    2.09KB        0
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
#> 1 sset(x, x %in_% y)   81.9µs    119µs     8606.    83.2KB     6.54
#> 2 sset(x, x %in% y)   156.4µs    238µs     4096.   285.4KB     8.90
#> 3 x[x %in% y]         131.2µs    210µs     4802.   324.5KB    13.8
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            137.5µs  527.1µs     1899.    1.53MB    27.5 
#> 2 df[0:10^5, , drop = FALSE]   6.01ms   7.14ms      141.    4.83MB     7.43
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
#> 1 sset(df, -10^4:0)             50.8ms   65.1ms     13.1      152MB     9.34
#> 2 df[-10^4:0, , drop = FALSE]  796.2ms  796.2ms      1.26     776MB     5.02
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
#> Warning: package 'data.table' was built under R version 4.4.1
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
#> 1 gcd(x)        1.3µs    1.4µs   658788.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       51.6ms   51.9ms      19.2        0B        0
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
x <- rep(TRUE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   2.52ms   3.92ms      252.    3.81MB     4.31
#> 2 base_which    654.2µs   2.68ms      379.    7.63MB    12.0
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    741µs    753µs     1311.        0B      0  
#> 2 base_which      453µs    464µs     2081.    3.81MB     28.7
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   2.02ms   2.71ms      364.    1.91MB     2.05
#> 2 base_which    776.9µs   1.77ms      558.    7.63MB    17.3
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   3.63ms   4.93ms      203.    3.81MB     2.07
#> 2 base_which    913.1µs   2.96ms      332.    9.54MB    11.8
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   2.42ms   3.01ms      331.    1.89MB     2.08
#> 2 base_which     3.16ms   4.09ms      245.     5.7MB     6.68
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
#> 1 cheapr_factor   9.19ms   9.91ms     97.6     4.59MB        0
#> 2 base_factor   529.85ms 529.85ms      1.89   27.84MB        0
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   4.34ms   4.93ms    202.      1.53MB     2.02
#> 2 base_factor   839.96ms 839.96ms      1.19   22.79MB     0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 221.21ms 224.83ms     4.44     5.23MB        0
#> 2 base_factor      3.09s    3.09s     0.324   54.35MB        0
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   4.88ms    6.5ms     152.     3.49MB     2.21
#> 2 base_factor    54.53ms   60.2ms      16.7   39.89MB     0
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
#> 1 cheapr_intersect   2.71ms   3.37ms      292.    1.18MB     2.17
#> 2 base_intersect     4.17ms   5.01ms      196.    5.16MB     2.18
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   2.99ms   3.12ms      312.    1.76MB     0   
#> 2 base_setdiff     4.45ms   5.43ms      185.    5.71MB     2.17
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.75ms   1.81ms      535.  781.34KB     2.17
#> 2 base          2.5ms   2.96ms      336.    2.53MB     2.18
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.75ms   1.83ms      533.  787.84KB     0   
#> 2 base         2.67ms   3.14ms      313.    2.91MB     2.17
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
#> 1 cheapr_cut    210ms    211ms      4.73    38.2MB     0   
#> 2 base_cut      487ms    509ms      1.96   267.1MB     2.95
```
