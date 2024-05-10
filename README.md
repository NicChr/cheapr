
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
#> 1 "na_locf(x)"                     984.2µs  992.5µs     1000.        0B      0  
#> 2 "vec_fill_missing(x, direction…   4.63ms   4.94ms      201.    11.4MB     67.1
mark(na_locf(x), vec_fill_missing(x, direction = "down"))
#> # A tibble: 2 × 6
#>   expression                           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 "na_locf(x)"                     986.3µs  990.1µs      997.        0B       0 
#> 2 "vec_fill_missing(x, direction…   3.85ms   5.19ms      190.    11.4MB     136.
```

All the `NA` handling functions in cheapr can make use of multiple cores
on your machine using openMP.

``` r
# 1 core by default
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)        982µs    986µs     1009.        0B      0  
#> 2 sum(is.na(x))    777µs      2ms      522.    3.81MB     52.7
# 4 cores
options(cheapr.cores = 4)
mark(num_na(x), sum(is.na(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 num_na(x)      256.3µs    298µs     3115.        0B      0  
#> 2 sum(is.na(x))    1.6ms      2ms      492.    3.81MB     44.7
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
#> 1 row_na_counts(m)    1.06ms   2.26ms      446.   12.88KB      0  
#> 2 rowSums(is.na(m))   2.61ms   3.85ms      264.    3.82MB     25.0
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)   679.7µs  790.3µs     1223.   12.88KB      0  
#> 2 colSums(is.na(m))   2.65ms   3.05ms      323.    3.82MB     32.5
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6)
x[sample.int(10^6, 10^5)] <- NA
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      777µs   1.97ms      502.    3.81MB     98.8
#> 2 is_na(x)      419µs  905.7µs     1101.    3.82MB     98.9

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
#> 1 is.na(hours)    1.19s    1.19s     0.840   61.05MB    0.840
#> 2 is_na(hours)   5.02ms   6.27ms   154.       9.97MB    9.97
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
#> 1   x integer         0          1      100 50.49     1    25  50   76  100
#> 2   z numeric         0          1 10000000     0 -5.55 -0.67   0 0.67 5.67
#>    iqr    sd  hist
#> 1   51 28.87 ▇▇▇▇▇
#> 2 1.35     1 ▁▂▇▁▁
#> 
#> ----- Categorical -----
#>   col  class n_missing p_complete n_unique n_levels min max
#> 1   y factor         0          1       26       26   A   Z
mark(overview(df))
#> # A tibble: 1 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 overview(df)    971ms    971ms      1.03    76.3MB        0
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
#> 1 sset(x, x %in_% y)     80µs    119µs     8722.    88.3KB     4.24
#> 2 sset(x, x %in% y)     163µs    238µs     4200.   285.4KB     4.26
#> 3 x[x %in% y]           130µs    208µs     4864.   324.5KB     9.04
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            235.8µs  403.7µs     2378.    1.53MB    16.1 
#> 2 df[0:10^5, , drop = FALSE]   6.11ms   7.46ms      130.    4.83MB     4.20
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
#> 1 sset(df, -10^4:0)             29.9ms   40.5ms     19.7      152MB    13.8 
#> 2 df[-10^4:0, , drop = FALSE]  732.9ms  732.9ms      1.36     776MB     6.82
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

lag_(1:10, 3) # Lag(3)
#>  [1] NA NA NA  1  2  3  4  5  6  7
lag_(1:10, -3) # Lead(3)
#>  [1]  4  5  6  7  8  9 10 NA NA NA

# Using an example from data.table
library(data.table)
dt <- data.table(year=2010:2014, v1=runif(5), v2=1:5, v3=letters[1:5])

# Similar to data.table::shift()

lag_(dt, 1) # Lag 
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:    NA        NA    NA   <NA>
#> 2:  2010 0.4443658     1      a
#> 3:  2011 0.9752701     2      b
#> 4:  2012 0.7064650     3      c
#> 5:  2013 0.8423412     4      d
lag_(dt, -1) # Lead
#>     year         v1    v2     v3
#>    <int>      <num> <int> <char>
#> 1:  2011 0.97527012     2      b
#> 2:  2012 0.70646498     3      c
#> 3:  2013 0.84234116     4      d
#> 4:  2014 0.08024799     5      e
#> 5:    NA         NA    NA   <NA>
```

With `lag_` we can update variables by reference, including entire data
frames

``` r
# At the moment, shift() cannot do this
lag_(dt, set = TRUE)
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:    NA        NA    NA   <NA>
#> 2:  2010 0.4443658     1      a
#> 3:  2011 0.9752701     2      b
#> 4:  2012 0.7064650     3      c
#> 5:  2013 0.8423412     4      d

dt # Was updated by reference
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:    NA        NA    NA   <NA>
#> 2:  2010 0.4443658     1      a
#> 3:  2011 0.9752701     2      b
#> 4:  2012 0.7064650     3      c
#> 5:  2013 0.8423412     4      d
```

`lag2_` is a more generalised variant that supports vectors of lags,
custom ordering and run lengths.

``` r
lag2_(dt, order = 5:1) # Reverse order lag (same as lead)
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:  2010 0.4443658     1      a
#> 2:  2011 0.9752701     2      b
#> 3:  2012 0.7064650     3      c
#> 4:  2013 0.8423412     4      d
#> 5:    NA        NA    NA   <NA>
lag2_(dt, -1) # Same as above
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:  2010 0.4443658     1      a
#> 2:  2011 0.9752701     2      b
#> 3:  2012 0.7064650     3      c
#> 4:  2013 0.8423412     4      d
#> 5:    NA        NA    NA   <NA>
lag2_(dt, c(1, -1)) # Alternating lead/lag
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:    NA        NA    NA   <NA>
#> 2:  2011 0.9752701     2      b
#> 3:  2010 0.4443658     1      a
#> 4:  2013 0.8423412     4      d
#> 5:  2012 0.7064650     3      c
lag2_(dt, c(-1, 0, 0, 0, 0)) # Lead e.g. only first row
#>     year        v1    v2     v3
#>    <int>     <num> <int> <char>
#> 1:  2010 0.4443658     1      a
#> 2:  2010 0.4443658     1      a
#> 3:  2011 0.9752701     2      b
#> 4:  2012 0.7064650     3      c
#> 5:  2013 0.8423412     4      d
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
#> 1 gcd(x)        1.3µs    1.4µs   655072.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       54.7ms   54.8ms      18.2        0B        0
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
#> 1 cheapr_which   2.06ms   3.65ms      277.    3.81MB     2.06
#> 2 base_which    698.8µs   2.61ms      387.    7.63MB     6.71
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    218µs    252µs     3637.        0B      0  
#> 2 base_which      454µs    462µs     2147.    3.81MB     17.1
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.26ms   2.06ms      477.    1.91MB     2.04
#> 2 base_which    785.3µs   1.77ms      567.    7.63MB    11.3
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   3.03ms    4.3ms      233.    3.81MB     2.08
#> 2 base_which      891µs   2.86ms      354.    9.54MB     6.56
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.88ms   2.68ms      380.    1.89MB     2.07
#> 2 base_which      3.2ms   4.06ms      250.    5.71MB     4.30
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
#> 1 cheapr_factor    9.2ms   9.95ms     94.2     4.59MB        0
#> 2 base_factor    611.6ms 611.61ms      1.64   27.84MB        0
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   4.36ms   5.29ms    187.      1.53MB        0
#> 2 base_factor   919.07ms 919.07ms      1.09   22.79MB        0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor 220.09ms 220.12ms     4.54     5.23MB    0    
#> 2 base_factor      3.17s    3.17s     0.316   54.35MB    0.316
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   5.27ms     11ms      97.4    3.49MB     0   
#> 2 base_factor    55.12ms   61.2ms      16.5   39.89MB     2.06
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
#> 1 cheapr_intersect    2.8ms   2.92ms      331.    1.18MB     0   
#> 2 base_intersect     4.45ms   5.24ms      182.    5.16MB     2.16
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   3.02ms   3.13ms      317.    1.76MB     0   
#> 2 base_setdiff     4.71ms   5.32ms      183.    5.71MB     2.15
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.75ms   1.81ms      544.  781.34KB     0   
#> 2 base         2.61ms   2.97ms      331.    2.53MB     2.15
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.71ms   1.86ms      517.  787.84KB     0   
#> 2 base         2.73ms   3.05ms      322.    2.91MB     2.15
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
#> 1 cheapr_cut    142ms    143ms      7.01    38.1MB     0   
#> 2 base_cut      480ms    499ms      2.00   267.1MB     3.00
```
