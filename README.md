
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

### Scalars and `NA`

Because R mostly uses vectors and vectorised operations, this means that
there are few scalar-optimised operations.

cheapr provides tools to efficiently count, find, replace and remove
scalars.

``` r
# Setup data with NA values
set.seed(42)
x <- sample(1:5, 30, TRUE)
x <- na_insert(x, n = 7)

counts(x) # Fast counts
#>   key count
#> 1   1     6
#> 2   5     4
#> 3  NA     7
#> 4   2     6
#> 5   4     4
#> 6   3     3
```

`NA` functions

``` r
na_count(x)
#> [1] 7
na_rm(x)
#>  [1] 1 5 1 2 4 2 1 4 5 4 2 3 1 1 3 4 5 5 2 3 2 1 2
na_find(x)
#> [1]  4  8 11 15 22 24 26
na_replace(x, -99)
#>  [1]   1   5   1 -99   2   4   2 -99   1   4 -99   5   4   2 -99   3   1   1   3
#> [20]   4   5 -99   5 -99   2 -99   3   2   1   2
```

Scalar functions

``` r
val_count(x, 3)
#> [1] 3
val_rm(x, 3)
#>  [1]  1  5  1 NA  2  4  2 NA  1  4 NA  5  4  2 NA  1  1  4  5 NA  5 NA  2 NA  2
#> [26]  1  2
val_find(x, 3)
#> [1] 16 19 27
val_replace(x, 3, 99)
#>  [1]  1  5  1 NA  2  4  2 NA  1  4 NA  5  4  2 NA 99  1  1 99  4  5 NA  5 NA  2
#> [26] NA 99  2  1  2
```

Scalar based case-match

``` r
val_match(
  x, 
  1 ~ "one", 
  2 ~ "two",
  3 ~ "three", 
  .default = ">3"
)
#>  [1] "one"   ">3"    "one"   ">3"    "two"   ">3"    "two"   ">3"    "one"  
#> [10] ">3"    ">3"    ">3"    ">3"    "two"   ">3"    "three" "one"   "one"  
#> [19] "three" ">3"    ">3"    ">3"    ">3"    ">3"    "two"   ">3"    "three"
#> [28] "two"   "one"   "two"
```

## Efficient NA counts by row/col

``` r
m <- matrix(na_insert(rnorm(10^6), prop = 1/4), ncol = 10^3)
# Number of NA values by row
mark(row_na_counts(m), 
     rowSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 row_na_counts(m)   456.6µs  488.7µs     1906.   13.14KB      0  
#> 2 rowSums(is.na(m))   2.78ms   3.07ms      303.    3.85MB     29.6
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)    1.35ms   1.46ms      665.   13.14KB      0  
#> 2 colSums(is.na(m))   1.26ms    1.5ms      609.    3.82MB     61.9
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6) |> 
  na_insert(10^5)
options(cheapr.cores = 4)
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      589µs    738µs     1320.    3.81MB     229.
#> 2 is_na(x)      146µs    173µs     5278.    3.82MB     530.
options(cheapr.cores = 1)
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      571µs    650µs     1460.    3.81MB     182.
#> 2 is_na(x)      376µs    446µs     2188.    3.81MB     263.

### posixlt method is much faster
hours <- as.POSIXlt(seq.int(0, length.out = 10^6, by = 3600),
                    tz = "UTC") |> 
  na_insert(10^5)

mark(is.na(hours), is_na(hours))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(hours)    1.04s    1.04s     0.963   61.05MB    0.963
#> 2 is_na(hours)   3.87ms   4.06ms   225.       7.65MB   23.9
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
df <- new_df(x = c(1, NA, 3),
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
df <- new_df(
  x = sample.int(100, 10^6, TRUE),
  y = as_factor(sample(LETTERS, 10^6, TRUE)),
  z = rnorm(10^6)
)
overview(df)
#> obs: 1000000 
#> cols: 3 
#> 
#> ----- Numeric -----
#>   col n_missng p_complt n_unique     mean    p0   p25      p50   p75   p100
#> 1   x        0        1      100    50.52     1    25       51    76    100
#> 2   z        0        1  1000000 -0.00038 -4.58 -0.67 -0.00062  0.68   5.08
#>     iqr    sd  hist
#> 1    51 28.88 ▇▇▇▇▇
#> 2  1.35     1 ▁▃▇▂▁
#> 
#> ----- Categorical -----
#>   col n_missng p_complt n_unique n_levels min max
#> 1   y        0        1       26       26   A   Z
mark(overview(df, hist = FALSE))
#> # A tibble: 1 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 overview(df, hist = FALSE)   69.7ms   71.1ms      13.9        0B        0
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
#> 1 sset(x, x %in_% y)   87.1µs    105µs     8617.      86KB     10.4
#> 2 sset(x, x %in% y)   147.7µs    162µs     5721.     286KB     25.8
#> 3 x[x %in% y]         142.5µs    157µs     5910.     325KB     32.2
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            122.8µs  188.5µs     4041.    1.53MB    69.1 
#> 2 df[0:10^5, , drop = FALSE]   6.15ms   6.91ms      141.    4.83MB     4.15
```

It also accepts negative indexes

``` r
mark(sset(df, -10^4:0), 
     df[-10^4:0, , drop = FALSE],
     check = FALSE) # The only difference is the row names
#> # A tibble: 2 × 6
#>   expression                       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, -10^4:0)            726.1µs   2.02ms     537.     15.1MB     131.
#> 2 df[-10^4:0, , drop = FALSE]   17.8ms  20.15ms      49.6    72.5MB     447.
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
#> Error in sset(x, TRUE): `length(i)` must match `length(x)` when `i` is a logical vector

# This is equivalent 
x[TRUE]
#> [1]  1  5 NA NA -5
# to..
sset(x)
#> [1]  1  5 NA NA -5
```

## Combining vectors fast and consistently

``` r
x <- as_factor(letters)
```

Base R combining

``` r
c(x, letters)
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
#> [16] "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "a"  "b"  "c"  "d" 
#> [31] "e"  "f"  "g"  "h"  "i"  "j"  "k"  "l"  "m"  "n"  "o"  "p"  "q"  "r"  "s" 
#> [46] "t"  "u"  "v"  "w"  "x"  "y"  "z"
c(letters, x)
#>  [1] "a"  "b"  "c"  "d"  "e"  "f"  "g"  "h"  "i"  "j"  "k"  "l"  "m"  "n"  "o" 
#> [16] "p"  "q"  "r"  "s"  "t"  "u"  "v"  "w"  "x"  "y"  "z"  "1"  "2"  "3"  "4" 
#> [31] "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19"
#> [46] "20" "21" "22" "23" "24" "25" "26"
```

With cheapr the order of arguments doesn’t affect the outcome type

``` r
c_(x, letters);c_(letters, x)
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z a b c d e f g h i j k l
#> [39] m n o p q r s t u v w x y z
#> Levels: a b c d e f g h i j k l m n o p q r s t u v w x y z
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z a b c d e f g h i j k l
#> [39] m n o p q r s t u v w x y z
#> Levels: a b c d e f g h i j k l m n o p q r s t u v w x y z
```

Same goes for other types likes Dates and Date-Times

``` r
today <- Sys.Date()
now <- Sys.time()

c(today, now);c(now, today) # base
#> [1] "2025-11-08" "2025-11-08"
#> [1] "2025-11-08 06:25:02 GMT" "2025-11-08 00:00:00 GMT"
c_(today, now);c_(now, today) # cheapr
#> [1] "2025-11-08 00:00:00 GMT" "2025-11-08 06:25:02 GMT"
#> [1] "2025-11-08 06:25:02 GMT" "2025-11-08 00:00:00 GMT"
```

`c_()` combines date frames by row

``` r
sset(iris, 1) |> 
  c_(sset(iris, 2))
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 1          5.1         3.5          1.4         0.2  setosa
#> 2          4.9         3.0          1.4         0.2  setosa
```

Alternatively combine by column

``` r
sset(iris, 1:3, j = 1) |> 
  col_c(sset(iris, 1:3, j = 2))
#>   Sepal.Length Sepal.Width
#> 1          5.1         3.5
#> 2          4.9         3.0
#> 3          4.7         3.2
```

## Casting and coercion

We can cast from one type to another with `cast()`

``` r
ints <- 1:10
dbls <- seq_(from = 1, to = 10, by = 0.5)
chrs <- letters
fctr <- as_factor(letters)
df <- new_df(a = ints, b = dbls, c = chrs, d = fctr)

cast(ints, dbls) |> print() |> class()
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> [1] "numeric"
cast(dbls, ints) |> print() |> class()
#>  [1]  1  1  2  2  3  3  4  4  5  5  6  6  7  7  8  8  9  9 10
#> [1] "integer"
cast(dbls, chrs) |> print() |> class()
#>  [1] "1"   "1.5" "2"   "2.5" "3"   "3.5" "4"   "4.5" "5"   "5.5" "6"   "6.5"
#> [13] "7"   "7.5" "8"   "8.5" "9"   "9.5" "10"
#> [1] "character"
cast(chrs, fctr) |> print() |> class()
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z
#> Levels: a b c d e f g h i j k l m n o p q r s t u v w x y z
#> [1] "factor"
cast(fctr, chrs) |> print() |> class()
#>  [1] "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "o" "p" "q" "r" "s"
#> [20] "t" "u" "v" "w" "x" "y" "z"
#> [1] "character"
cast(dbls, fctr) |> print() |> class()
#>  [1] <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> [16] <NA> <NA> <NA> <NA>
#> Levels: a b c d e f g h i j k l m n o p q r s t u v w x y z
#> [1] "factor"
cast(fctr, dbls) |> print() |> class()
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> [26] 26
#> [1] "numeric"
cast(ints, df) |> print() |> class()
#>    value
#> 1      1
#> 2      2
#> 3      3
#> 4      4
#> 5      5
#> 6      6
#> 7      7
#> 8      8
#> 9      9
#> 10    10
#> [1] "data.frame"
```

We can also cast multiple objects to a common type

``` r
cast_common(ints, dbls)
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[2]]
#>  [1]  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0  5.5  6.0  6.5  7.0  7.5  8.0
#> [16]  8.5  9.0  9.5 10.0
cast_common(ints, dbls, chrs)
#> [[1]]
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"
#> 
#> [[2]]
#>  [1] "1"   "1.5" "2"   "2.5" "3"   "3.5" "4"   "4.5" "5"   "5.5" "6"   "6.5"
#> [13] "7"   "7.5" "8"   "8.5" "9"   "9.5" "10" 
#> 
#> [[3]]
#>  [1] "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "o" "p" "q" "r" "s"
#> [20] "t" "u" "v" "w" "x" "y" "z"
cast_common(ints, dbls, chrs, fctr)
#> [[1]]
#>  [1] 1  2  3  4  5  6  7  8  9  10
#> 45 Levels: 1 2 3 4 5 6 7 8 9 10 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 a b c ... z
#> 
#> [[2]]
#>  [1] 1   1.5 2   2.5 3   3.5 4   4.5 5   5.5 6   6.5 7   7.5 8   8.5 9   9.5 10 
#> 45 Levels: 1 2 3 4 5 6 7 8 9 10 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 a b c ... z
#> 
#> [[3]]
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z
#> 45 Levels: 1 2 3 4 5 6 7 8 9 10 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 a b c ... z
#> 
#> [[4]]
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z
#> 45 Levels: 1 2 3 4 5 6 7 8 9 10 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 a b c ... z
cast_common(ints, dbls, chrs, fctr, df)
#> [[1]]
#>    value
#> 1      1
#> 2      2
#> 3      3
#> 4      4
#> 5      5
#> 6      6
#> 7      7
#> 8      8
#> 9      9
#> 10    10
#> 
#> [[2]]
#>    value
#> 1    1.0
#> 2    1.5
#> 3    2.0
#> 4    2.5
#> 5    3.0
#> 6    3.5
#> 7    4.0
#> 8    4.5
#> 9    5.0
#> 10   5.5
#> 11   6.0
#> 12   6.5
#> 13   7.0
#> 14   7.5
#> 15   8.0
#> 16   8.5
#> 17   9.0
#> 18   9.5
#> 19  10.0
#> 
#> [[3]]
#>    value
#> 1      a
#> 2      b
#> 3      c
#> 4      d
#> 5      e
#> 6      f
#> 7      g
#> 8      h
#> 9      i
#> 10     j
#> 11     k
#> 12     l
#> 13     m
#> 14     n
#> 15     o
#> 16     p
#> 17     q
#> 18     r
#> 19     s
#> 20     t
#> 21     u
#> 22     v
#> 23     w
#> 24     x
#> 25     y
#> 26     z
#> 
#> [[4]]
#>    value
#> 1      a
#> 2      b
#> 3      c
#> 4      d
#> 5      e
#> 6      f
#> 7      g
#> 8      h
#> 9      i
#> 10     j
#> 11     k
#> 12     l
#> 13     m
#> 14     n
#> 15     o
#> 16     p
#> 17     q
#> 18     r
#> 19     s
#> 20     t
#> 21     u
#> 22     v
#> 23     w
#> 24     x
#> 25     y
#> 26     z
#> 
#> [[5]]
#>     a    b c d
#> 1   1  1.0 a a
#> 2   2  1.5 b b
#> 3   3  2.0 c c
#> 4   4  2.5 d d
#> 5   5  3.0 e e
#> 6   6  3.5 f f
#> 7   7  4.0 g g
#> 8   8  4.5 h h
#> 9   9  5.0 i i
#> 10 10  5.5 j j
#> 11  1  6.0 k k
#> 12  2  6.5 l l
#> 13  3  7.0 m m
#> 14  4  7.5 n n
#> 15  5  8.0 o o
#> 16  6  8.5 p p
#> 17  7  9.0 q q
#> 18  8  9.5 r r
#> 19  9 10.0 s s
#> 20 10  1.0 t t
#> 21  1  1.5 u u
#> 22  2  2.0 v v
#> 23  3  2.5 w w
#> 24  4  3.0 x x
#> 25  5  3.5 y y
#> 26  6  4.0 z z
```

When common-casting factors, their levels are combined

``` r
cast_common(fctr, as_factor(LETTERS))
#> [[1]]
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z
#> 52 Levels: a b c d e f g h i j k l m n o p q r s t u v w x y z A B C D E ... Z
#> 
#> [[2]]
#>  [1] A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
#> 52 Levels: a b c d e f g h i j k l m n o p q r s t u v w x y z A B C D E ... Z
cast_common(as_factor(LETTERS), fctr)
#> [[1]]
#>  [1] A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
#> 52 Levels: A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e ... z
#> 
#> [[2]]
#>  [1] a b c d e f g h i j k l m n o p q r s t u v w x y z
#> 52 Levels: A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e ... z
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
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:cheapr':
#> 
#>     address
data.table::setDTthreads(1);
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
#> 1 gcd(x)        800ns    900ns   992704.        0B        0
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       30.7ms   32.1ms      31.2        0B        0
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

Use `as_list` to return a list of sequences

``` r
seq_(start, end, increments, as_list = TRUE)
#> [[1]]
#> [1] 1 2 3 4 5
#> 
#> [[2]]
#> [1] 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
#> 
#> [[3]]
#>  [1] 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8
#> [20] 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7
#> [39] 4.8 4.9 5.0
```

If you know the sizes of your sequences beforehand, use `sequence_()`

``` r
seq_sizes <- c(3, 5, 10)
sequence_(seq_sizes, from = 0, by = 1/3, as_list = TRUE)
#> [[1]]
#> [1] 0.0000000 0.3333333 0.6666667
#> 
#> [[2]]
#> [1] 0.0000000 0.3333333 0.6666667 1.0000000 1.3333333
#> 
#> [[3]]
#>  [1] 0.0000000 0.3333333 0.6666667 1.0000000 1.3333333 1.6666667 2.0000000
#>  [8] 2.3333333 2.6666667 3.0000000
```

You can also calculate sequence sizes, starts, ends and increments

``` r
seq_size(from = 1, to = 10, by = c(0.5, 1))
#> [1] 19 10
seq_start(size = c(19, 10), to = 10, by = c(0.5, 1))
#> [1] 1 1
seq_end(size = c(19, 10), from = 1, by = c(0.5, 1))
#> [1] 10 10
seq_increment(size = c(19, 10), from = 1, to = 10)
#> [1] 0.5 1.0
```

## Math in-place

cheapr provides a full set of common math functions that can transform
numeric vectors in-place (no copies)

``` r
(x <- seq(0, 5, by = 0.5))
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0

# x is modified in-place
set_add(x, 10);x
#>  [1] 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5 15.0
#>  [1] 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5 15.0
set_subtract(x, 10);x
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
set_multiply(x, 10);x
#>  [1]  0  5 10 15 20 25 30 35 40 45 50
#>  [1]  0  5 10 15 20 25 30 35 40 45 50
set_divide(x, 10);x
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0

set_change_sign(x);x
#>  [1]  0.0 -0.5 -1.0 -1.5 -2.0 -2.5 -3.0 -3.5 -4.0 -4.5 -5.0
#>  [1]  0.0 -0.5 -1.0 -1.5 -2.0 -2.5 -3.0 -3.5 -4.0 -4.5 -5.0
set_abs(x);x
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
set_round(x);x
#>  [1] 0 0 1 2 2 2 3 4 4 4 5
#>  [1] 0 0 1 2 2 2 3 4 4 4 5
set_log(x);x
#>  [1]      -Inf      -Inf 0.0000000 0.6931472 0.6931472 0.6931472 1.0986123
#>  [8] 1.3862944 1.3862944 1.3862944 1.6094379
#>  [1]      -Inf      -Inf 0.0000000 0.6931472 0.6931472 0.6931472 1.0986123
#>  [8] 1.3862944 1.3862944 1.3862944 1.6094379
```

These in-place functions are not always faster than using normal R math
functions. This becomes apparent when performing multiple operations
which R can process simultaneously.

``` r
x <- rnorm(10^6)
mark(
  x * 10 * 20 + 1 - 1 , 
  set_subtract(set_add(set_multiply(set_multiply(x, 10), 20), 1), 1)
)
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 x * 10 * 20 + 1 - 1                 1.37ms 1.75ms      503.    7.63MB     43.4
#> 2 set_subtract(set_add(set_multiply(… 3.14ms 3.34ms      287.        0B      0
```

### `.args`

cheapr now provides `.args` as a means of providing a list of arguments
instead of `...`. This is designed to replace the use of `do.call()`.

In practice this means that users can either supply objects directly to
the dots `...` or as a list of objects.

``` r
# The below lines are equivalent
c_(1, 2, 3)
#> [1] 1 2 3
c_(.args = list(1, 2, 3))
#> [1] 1 2 3
```

A very common scenario is having a list of objects that you would like
to combine into a vector. Normally one would call `do.call(c, x)` but it
is much more efficient to use the `.args` argument in `c_()`.

``` r
x <- rep(list(0), 10^5)

mark(
  do.call(c, x),
  c_(.args = x)
)
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 do.call(c, x)   2.73ms   3.52ms      268.     781KB   147.  
#> 2 c_(.args = x)   1.04ms   1.15ms      844.     781KB     4.16

# Matches the speed of `unlist()` without removing attributes
unlist(list(Sys.Date()), recursive = FALSE)
#> [1] 20400
c_(.args = list(Sys.Date()))
#> [1] "2025-11-08"
```

## Recycling

Fast base-R style recycling using `recycle()`

``` r
recycle(letters, pi)
#> [[1]]
#>  [1] "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "o" "p" "q" "r" "s"
#> [20] "t" "u" "v" "w" "x" "y" "z"
#> 
#> [[2]]
#>  [1] 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593
#>  [9] 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593
#> [17] 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593
#> [25] 3.141593 3.141593

# Data frame rows are recycled
recycle(vector = 1:10, data = cars)
#> $vector
#>  [1]  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5
#> [26]  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10
#> 
#> $data
#>    speed dist
#> 1      4    2
#> 2      4   10
#> 3      7    4
#> 4      7   22
#> 5      8   16
#> 6      9   10
#> 7     10   18
#> 8     10   26
#> 9     10   34
#> 10    11   17
#> 11    11   28
#> 12    12   14
#> 13    12   20
#> 14    12   24
#> 15    12   28
#> 16    13   26
#> 17    13   34
#> 18    13   34
#> 19    13   46
#> 20    14   26
#> 21    14   36
#> 22    14   60
#> 23    14   80
#> 24    15   20
#> 25    15   26
#> 26    15   54
#> 27    16   32
#> 28    16   40
#> 29    17   32
#> 30    17   40
#> 31    17   50
#> 32    18   42
#> 33    18   56
#> 34    18   76
#> 35    18   84
#> 36    19   36
#> 37    19   46
#> 38    19   68
#> 39    20   32
#> 40    20   48
#> 41    20   52
#> 42    20   56
#> 43    20   64
#> 44    22   66
#> 45    23   54
#> 46    24   70
#> 47    24   92
#> 48    24   93
#> 49    24  120
#> 50    25   85

# Using .args
recycle(.args = list(letters, pi))
#> [[1]]
#>  [1] "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "o" "p" "q" "r" "s"
#> [20] "t" "u" "v" "w" "x" "y" "z"
#> 
#> [[2]]
#>  [1] 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593
#>  [9] 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593
#> [17] 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593 3.141593
#> [25] 3.141593 3.141593
```

Sizes are recycled to the common maximum, except when a vector is length
0 (excluding NULL which is ignored), in which case they are all recycled
to length 0.

``` r
recycle(a = 1:3, b = 1:10, c = iris, d = numeric())
#> $a
#> integer(0)
#> 
#> $b
#> integer(0)
#> 
#> $c
#> [1] Sepal.Length Sepal.Width  Petal.Length Petal.Width  Species     
#> <0 rows> (or 0-length row.names)
#> 
#> $d
#> numeric(0)
```

## Copying

cheapr provides some helpers in the form of `shallow_copy`, `semi_copy`
and `deep_copy`.

``` r
mark(shallow_copy(iris))
#> # A tibble: 1 × 6
#>   expression              min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>         <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 shallow_copy(iris)    300ns    400ns  1958902.    6.34KB        0
mark(deep_copy(iris))
#> # A tibble: 1 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 deep_copy(iris)    700ns    1.4µs   526576.    9.36KB     52.7
mark(semi_copy(iris))
#> # A tibble: 1 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 semi_copy(iris)    700ns    1.5µs   529008.    9.38KB        0
```

### `shallow_copy`

Shallow-copies list elements and attributes. When given an atomic vector
it full copies the vector and so is mostly useful for lists.

### `deep_copy`

Full (deep) copies everything, including attributes.

### `semi_copy`

Like `deep_copy` it deep-copies everything, excluding attributes, which
it shallow copies. In practice this turns out to be more efficient.

`semi_copy()` vs `deep_copy()`

``` r
df <- new_df(x = integer(10^6))
attr(df, "my_attr") <- integer(10^6)

# Take note of the memory allocation
mark(
  semi_copy(df), # Only deep copies the data
  deep_copy(df) # Deep copies "my_attr" as well
)
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 semi_copy(df)    114µs    128µs     6167.    3.81MB     266.
#> 2 deep_copy(df)    234µs    322µs     1960.    7.63MB     196.
```

### Attributes

With cheapr you can add and remove attributes flexibly using
`attrs_add()`.

To remove all attributes, use `attrs_rm()`.

To remove specific attributes, use `attrs_add(attr = NULL)`.

``` r
(x <- attrs_add(1:10, .length = 10, .type = "integer"))
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> attr(,".length")
#> [1] 10
#> attr(,".type")
#> [1] "integer"

attrs_add(x, .type = NULL) # Remove specific attribute '.type'
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> attr(,".length")
#> [1] 10
attrs_rm(x) # Clear all attributes
#>  [1]  1  2  3  4  5  6  7  8  9 10

# With .args
y <- 11:20
attrs_add(y, .args = attributes(x))
#>  [1] 11 12 13 14 15 16 17 18 19 20
#> attr(,".length")
#> [1] 10
#> attr(,".type")
#> [1] "integer"
```

Both functions allow setting attributes in-place. This turns out to be
very useful in avoiding implicit copies that R performs when it detects
that the data has been modified.

This must be used with care to not overwrite an existing object’s
attributes. Therefore it is best-practice to only use in-place attribute
manipulation on fresh objects, i.e objects that you can ensure are newly
created.

``` r
add_length_class <- function(x){
  attr(x, ".length") <- length(x)
  attr(x, ".class") <- class(x)
  x
}
add_length_class_in_place <- function(x){
  attrs_add(
    x, .length = length(x), .class = class(x),
    .set = TRUE
  )
}

# Notice the memory allocations
# we expect only 3.81 MB to be allocated
mark(
  add_length_class(integer(10^6)),
  add_length_class_in_place(integer(10^6)),
  iterations = 1
)
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                        <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 add_length_class(integer(10^6))     2.4ms   2.4ms      418.    3.81MB     418.
#> 2 add_length_class_in_place(intege… 364.1µs 364.1µs     2746.    3.81MB       0
mark(
  add_length_class(integer(10^6)),
  add_length_class_in_place(integer(10^6)),
  iterations = 1
)
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 add_length_class(integer(10^6))      322µs  322µs     3109.    3.81MB        0
#> 2 add_length_class_in_place(integer(1… 334µs  334µs     2990.    3.81MB        0
  

# R detected that the vector we created had been modified (because it was)
# and created a copy
# When we add the attributes in-place to our fresh object, no copies are
# made
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
#> 1 cheapr_which   2.55ms   3.13ms      293.    3.82MB     8.56
#> 2 base_which    559.1µs 654.45µs     1076.    7.63MB    88.7
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    118µs    124µs     7584.        0B       0 
#> 2 base_which      225µs    232µs     3998.    3.81MB     141.
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.34ms   1.56ms      572.    1.91MB     11.1
#> 2 base_which    490.5µs    752µs     1274.    7.63MB     86.9
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    895µs   1.15ms      873.    3.81MB     28.5
#> 2 base_which      724µs 887.05µs      839.    9.54MB     80.4
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which  577.4µs  647.4µs     1405.    1.89MB     24.8
#> 2 base_which     3.66ms   3.94ms      247.     5.7MB     13.1
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
#> 1 cheapr_factor   8.74ms   9.44ms    102.      6.13MB     7.14
#> 2 base_factor   300.12ms 300.12ms      3.33   27.84MB     3.33
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   2.98ms    3.2ms    296.      1.55MB     6.84
#> 2 base_factor    470.5ms  473.3ms      2.11   22.79MB     0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor  65.96ms  66.24ms    14.9      17.4MB     5.95
#> 2 base_factor      2.65s    2.65s     0.377    54.4MB     0
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   4.25ms   4.89ms     197.     3.49MB     9.49
#> 2 base_factor    43.35ms   46.5ms      21.8   39.89MB    36.3
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
#> 1 cheapr_intersect   2.05ms    2.2ms      421.    1.55MB     6.79
#> 2 base_intersect     3.93ms   4.52ms      218.    6.41MB    16.1
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   1.92ms   2.25ms      425.    2.15MB     9.50
#> 2 base_setdiff     3.96ms   4.34ms      221.    6.96MB    15.4
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.24ms   1.34ms      712.  781.34KB     4.28
#> 2 base          2.2ms   2.36ms      411.    2.53MB    11.9
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr        1.2ms   1.31ms      731.  792.32KB     6.76
#> 2 base         2.16ms   2.41ms      398.    2.91MB     9.20
```

### `as_discrete`

`as_discrete` is a cheaper alternative to `cut`

``` r
x <- rnorm(10^6)
b <- seq(0, max(x), 0.2)
mark(cheapr_cut = as_discrete(x, b, left = FALSE), 
     base_cut = cut(x, b))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_cut   14.2ms   14.8ms      66.6    3.92MB     4.44
#> 2 base_cut     25.2ms   25.5ms      38.2   15.32MB     7.63
```

### `if_else_`

A cheap alternative to `ifelse`

``` r
mark(
  if_else_(x >= 0, 1, -1),
  ifelse(x >= 0, 1, -1),
  data.table::fifelse(x >= 0, 1, -1)
)
#> # A tibble: 3 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                        <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 if_else_(x >= 0, 1, -1)            2.65ms  2.75ms     348.     11.4MB     85.8
#> 2 ifelse(x >= 0, 1, -1)             17.04ms 17.92ms      55.4    53.4MB    200. 
#> 3 data.table::fifelse(x >= 0, 1, -…  5.36ms   5.6ms     172.     11.4MB     43.0
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
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 "case(x >= 0 ~ \"pos\", x < 0 ~ \"… 17.1ms 18.6ms      53.6    28.8MB     37.1
#> 2 "data.table::fcase(x >= 0, \"pos\"… 15.3ms 17.3ms      58.0    26.7MB     36.3
```

`val_match` is an even cheaper special variant of `case` when all LHS
expressions are length-1 vectors, i.e scalars

``` r
x <- round(rnorm(10^6))

mark(
  val_match(x, 1 ~ Inf, 2 ~ -Inf, .default = NaN),
     case(x == 1 ~ Inf, 
          x == 2 ~ -Inf, 
          .default = NaN),
     data.table::fcase(x == 1, Inf, 
          x == 2, -Inf, 
          rep_len(TRUE, length(x)), NaN)
     )
#> # A tibble: 3 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                        <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 val_match(x, 1 ~ Inf, 2 ~ -Inf, …  3.23ms  3.52ms     269.     8.79MB     36.9
#> 2 case(x == 1 ~ Inf, x == 2 ~ -Inf…  12.9ms 13.41ms      73.3   27.63MB     44.0
#> 3 data.table::fcase(x == 1, Inf, x…    11ms  12.1ms      81.2   30.52MB     81.2
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
#> 1 get_breaks(x, 20)   61.1µs   64.9µs    14803.        0B      0  
#> 2 pretty(x, 20)      411.3µs  538.1µs     1779.    1.91MB     36.5

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

<img src="man/figures/README-unnamed-chunk-52-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-52-2.png" width="100%" />
