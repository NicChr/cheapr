
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

cheapr is primarily an extension to R for developers that want to write
clean, fast and safe code without sacrificing any one of these core
design principles.

cheapr includes both an R API in the usual form of an R package, as well
as a C/C++ API for writing C/C++ code. It can be used interchangeably
with the R C API. C++17 or later is required to install and use cheapr.

**Important note:** While cheapr makes heavy usage of C++ templates and
variadic functions, it still internally uses the R C API only and throws
R errors via `Rf_error()`. This means that even though C++ templates are
heavily used, it still remains effectively a ‘pure C’ extension to the R
C API.

### Using C++ with cheapr

cheapr tries to strike a nice balance between speed, readability, safety
and flexibility. To achieve that speed and efficiency, we have opted to
internally rely on the R C API entirely.  
This doesn’t come without its tradeoffs, the biggest one being that
using cheapr with C++ or cpp11 becomes more complicated and requires
using `cheapr::r_safe`

### Notes on using cheapr with the R C API vs with C++/cpp11

<table style="width:92%;">
<colgroup>
<col style="width: 30%" />
<col style="width: 30%" />
<col style="width: 30%" />
</colgroup>
<thead>
<tr>
<th>Feature</th>
<th>R C API</th>
<th>C++/cpp11</th>
</tr>
</thead>
<tbody>
<tr>
<td>Protecting objects from GC</td>
<td><p>Pros: Fast</p>
<p>Cons: Manually protect with PROTECT/UNPROTECT</p></td>
<td>Pros: No manual protection needed.<br />
Cons: Slow if re-protecting an object many times (e.g. &gt;10k)</td>
</tr>
<tr>
<td>Safety</td>
<td>No need to use <code>r_safe</code> when not using C++ objects.</td>
<td>Must use <code>r_safe</code> on either cheapr functions or R C API
functions when C++ objects are created.</td>
</tr>
<tr>
<td>Type declarations</td>
<td>All R objects are declared as SEXP<br />
Pros: Sometimes this ambiguity is useful.<br />
Cons: Unclear what the return type should be.</td>
<td>Pros: Uses proper C++ types, e.g. integers, doubles, etc, making the
return type clear.<br />
Cons: Not flexible if needing an input which can be either integer or
double in which case cpp11::sexp can be used.</td>
</tr>
<tr>
<td>Exporting C/C++ functions to R</td>
<td>Use [[cpp11::register]]</td>
<td>Use [[cpp11::register]]</td>
</tr>
<tr>
<td>Using r_safe/cpp11::safe</td>
<td>No need to use <code>r_safe</code> when <strong>ONLY</strong> using
cheapr + R C API</td>
<td>You must use <code>r_safe</code> <strong>on</strong> cheapr/R C API
functions <strong>IF</strong> you have used any cpp11/C++ code.</td>
</tr>
</tbody>
</table>

When using C++ objects like for example std::string or cpp11::integers
you must call any pure R C function or any cheapr function with either
`cpp11::safe` or `cheapr::r_safe`. They both attempt to do the same
thing, catch any errors thrown by R and run any C++ destructors to make
sure no memory allocated to a C++ object is leaked. For more info see
cpp11’s FAQ section 15 [Should I call cpp11::unwind_protect()
manually?](https://cpp11.r-lib.org/articles/FAQ.html) The only
difference is that `r_safe` can accept template functions, which cheapr
makes heavy use of.

### DONT

``` cpp
#include <cpp11.hpp>
#include <cheapr_api.h>

cpp11::integers foo(cpp11::integers x){
  return cheapr::sset(x, 1, true); // First value of x
}
```

### DO

``` cpp
#include <cpp11.hpp>
#include <cheapr_api.h>

cpp11::integers foo(cpp11::integers x){
  return cheapr::r_safe(cheapr::sset)(x, 1, true); // First value of x
}
```

### OR (using C objects only)

``` cpp
#include <cheapr_api.h>

SEXP foo(SEXP x){
  return cheapr::sset(x, 1, true); // First value of x
}
```

### DONT

``` cpp
#include <cheapr_api.h>
SEXP foo(SEXP x){
  SEXP first = cheapr::sset(x, 1, true);
  first = coerce_vec(first, INTSXP);
  return first;
}
```

### DO

``` cpp
#include <cheapr_api.h>
SEXP foo(SEXP x){
  SEXP first = SHIELD(cheapr::sset(x, 1, true));
  first = SHIELD(coerce_vec(first, INTSXP));
  YIELD(2);
  return first;
}
```

## Exporting your C/C++ code to R

It is recommended to use the `[[cpp11::register]]` tag always to safely
export your C/C++ functions to R.

``` cpp
#include <cheapr_api.h>
[[cpp11::register]]
SEXP my_last(SEXP x){
  // max(length(x), 1)
  SEXP last_loc = SHIELD(as_r_vec(
    std::max(Rf_xlength(x), static_cast<R_xlen_t>(1))
  ));
  SEXP last_value = SHIELD(cheapr::sset(x, last_loc, true));
  YIELD(2);
  return last_value;
}
```

For help on getting started with C++ in R, see [Getting started with
cpp11](https://cpp11.r-lib.org/articles/cpp11.html)

### Using the cheapr C++ API

Let’s first load the required packages

``` r
library(cheapr) 
library(bench)
```

All the public user-facing C/C++ code is included in inst/include. To
make use of the API, simply include the cheapr API header file.

After this you have to link to cheapr either via the description file if
writing an R package

`LinkingTo: cheapr`

or by including the cpp11 tag `[[cpp11::linking_to("cheapr")]]` in your
C/C++ code.

``` r

setup_code <- '
#include <cheapr_api.h>
[[cpp11::linking_to("cheapr")]]
using namespace cheapr;

'
```

The functions can be found in the cheapr namespace

``` r
cpp11::cpp_source(
  code = paste_(
    setup_code,
    '
  [[cpp11::register]]
  bool foo(){
  return cheapr::is_r_na(NA_INTEGER);
  }
  '
  )
  , 
  cxx_std = "CXX17"
)
foo()
#> [1] TRUE
```

Write `using namespace cheapr` to make cheapr C++ fns available without
needing to use `cheapr::`

``` r
cpp11::cpp_source(
  code = paste_(
    setup_code,
    '
  using namespace cheapr;
  
  [[cpp11::register]]
  bool bar(){
  return is_r_na(NA_INTEGER);
  }
  '
  )
  , 
  cxx_std = "CXX17"
)
bar()
#> [1] TRUE
```

cheapr has many useful C++ functions you can use in your own C++ code.

One of the most powerful functions in cheapr’s C++ API is `new_r_vec`
which allows for creating R vectors by combining both C++ types as well
as other R vectors into one single R vector.

``` r
cpp11::cpp_source(
  code = paste_(
    setup_code,
    '
  [[cpp11::register]]
  SEXP foo(){
  return new_r_vec(1, 2, 3);
  }
  [[cpp11::register]]
  SEXP bar(){
  return new_r_vec(arg("x") = 1, arg("y") = 2, arg("z") = 3);
  }
  [[cpp11::register]]
  SEXP foobar(){
  return new_r_vec(1, R_NilValue, NA_INTEGER, "hello");
  }
  '
  )
  , 
  cxx_std = "CXX17"
)
foo()
#> [1] 1 2 3
bar()
#> x y z 
#> 1 2 3
foobar() # Implicit cast to character
#> [1] "1"     NA      "hello"
```

It acts like a C++ version of `base::c()` and utilises type-stable
common-casting. For more info on casting see the ‘Casting and Coercion’
section further below.

``` r
cpp11::cpp_source(
  code = paste_(
    setup_code,
    '
  [[cpp11::register]]
  SEXP foo(SEXP x, SEXP y){
  return new_r_vec(arg("x") = x, arg("y") = y);
  }
  '
  )
  , 
  cxx_std = "CXX17"
)
foo(pi, 1:3)
#>        x        y        y        y 
#> 3.141593 1.000000 2.000000 3.000000
```

Similarly, `new_r_list` can create an R list from both C++ types and R
vectors

``` r
cpp11::cpp_source(
  code = paste_(
    setup_code,
    '
  [[cpp11::register]]
  SEXP foo(){
  return new_r_list(
  true, // bool
  1, // int
  2.5, // Double
  3U, // Unsigned int
  4LL, // long-long int
  "hi", // const char *
  make_utf8_char("hey"), // CHARSXP
  make_utf8_str("hello") // STRSXP 
  );
  }
  '
  )
  , 
  cxx_std = "CXX17"
)
foo()
#> [[1]]
#> [1] TRUE
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 2.5
#> 
#> [[4]]
#> [1] 3
#> 
#> [[5]]
#> [1] 4
#> 
#> [[6]]
#> [1] "hi"
#> 
#> [[7]]
#> [1] "hey"
#> 
#> [[8]]
#> [1] "hello"
foo() |> c_(.args = _)
#> [1] "TRUE"  "1"     "2.5"   "3"     "4"     "hi"    "hey"   "hello"
```

Subsetting vectors with `sset()`

``` r
cpp11::cpp_source(
    code = paste_(
        setup_code,
        '
  [[cpp11::register]]
  SEXP foobar(SEXP x, SEXP i){
  return cheapr::sset(x, i, true);
  }
  '
    )
    , 
    cxx_std = "CXX17"
)
x <- 1:10
names(x) <- letters[1:10]
foobar(x, 3:1) # subset elements 3 to 1
#> c b a 
#> 3 2 1
foobar(x, "e") # Element with name "e"
#> e 
#> 5
foobar(x, -5)  # All elements except element 5
#>  a  b  c  d  f  g  h  i  j 
#>  1  2  3  4  6  7  8  9 10
foobar(x, c(0, NA_integer_, 100)) # Elements that don't exist return NA
#> <NA> <NA> 
#>   NA   NA
```

Repeating vectors with `rep_len()`, `rep_()` and `rep_each()`

``` r
cpp11::cpp_source(
    code = paste_(
        setup_code,
        '
  [[cpp11::register]]
  SEXP cpp_rep(SEXP x, SEXP times){
  return cheapr::rep(x, times);
  }
    [[cpp11::register]]
  SEXP cpp_rep_len(SEXP x, int64_t n){
  return cheapr::rep_len(x, n);
  }
    [[cpp11::register]]
  SEXP cpp_rep_each(SEXP x, SEXP each){
  return cheapr::rep_each(x, each);
  }
  '
    )
    , 
    cxx_std = "CXX17"
)
x <- 1:10
cpp_rep(x, 3)
#>  [1]  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5
#> [26]  6  7  8  9 10
cpp_rep_len(x, 20)
#>  [1]  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10
cpp_rep_each(x, 3)
#>  [1]  1  1  1  2  2  2  3  3  3  4  4  4  5  5  5  6  6  6  7  7  7  8  8  8  9
#> [26]  9  9 10 10 10
```

There are many more useful C++ functions in the API. Navigate to
inst/include to see all of them.

### cheapr R API

Some common R operations that cheapr can do much faster and more
efficiently include:

- Handling `NA` values very efficiently

- Counting, finding, removing and replacing scalar values

- Type-stable one-way casting and common-casting

- Combining vectors

- Creating and manipulating factors

- Pasting strings

- Creating multiple sequences in a vectorised way

- Sub-setting vectors and data frames efficiently

- Safe, flexible and fast greatest common divisor and lowest common
  multiple

- Lags/leads

- `integer64` support

- In-memory Math (no copies, vectors updated by reference)

- Summary statistics

- Counts

- Modifying lists

- Recycling

- Binning of continuous data

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
#> 1 row_na_counts(m)   455.7µs  468.4µs     1691.   13.09KB      0  
#> 2 rowSums(is.na(m))    2.2ms   3.49ms      263.    3.85MB     27.1
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)    1.36ms   1.53ms      558.   13.09KB      0  
#> 2 colSums(is.na(m))   1.27ms   1.46ms      600.    3.82MB     60.7
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
#> 1 is.na(x)      568µs    710µs     1244.    3.81MB     210.
#> 2 is_na(x)      159µs    193µs     4667.    3.82MB     465.
options(cheapr.cores = 1)
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      571µs    704µs     1300.    3.81MB     140.
#> 2 is_na(x)      358µs    422µs     2194.    3.81MB     220.

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
#> 1 is.na(hours)    1.14s    1.14s     0.873   61.05MB    0.873
#> 2 is_na(hours)   3.76ms   4.86ms   192.       7.65MB   20.0
```

It differs in 2 regards:

- List elements are regarded as `NA` only when that element is `NULL`
- For data frames, `is_na` returns a logical vector where `TRUE` defines
  an empty row of only `NA` values.

``` r
# List example
is.na(list(NA, NULL, 10))
#> [1]  TRUE FALSE FALSE
is_na(list(NA, NULL, 10))
#> [1] FALSE  TRUE FALSE

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
#> 1 overview(df, hist = FALSE)     79ms   84.4ms      11.8      512B        0
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
#> 1 sset(x, x %in_% y)   92.4µs    138µs     7080.      86KB     8.42
#> 2 sset(x, x %in% y)   145.7µs    176µs     4976.     286KB    21.5 
#> 3 x[x %in% y]         161.9µs    252µs     3687.     325KB    19.7
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            133.5µs  190.7µs     4228.    1.53MB   109.  
#> 2 df[0:10^5, , drop = FALSE]   6.41ms   7.56ms      123.    4.83MB     4.25
```

It also accepts negative indexes

``` r
mark(sset(df, -10^4:0), 
     df[-10^4:0, , drop = FALSE],
     check = FALSE) # The only difference is the row names
#> # A tibble: 2 × 6
#>   expression                       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, -10^4:0)              795µs    2.5ms     461.     15.1MB     131.
#> 2 df[-10^4:0, , drop = FALSE]   20.9ms     27ms      39.5    72.5MB     237.
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
#> [1] "2025-11-28" "2025-11-28"
#> [1] "2025-11-28 09:15:41 GMT" "2025-11-28 00:00:00 GMT"
c_(today, now);c_(now, today) # cheapr
#> [1] "2025-11-28 00:00:00 GMT" "2025-11-28 09:15:41 GMT"
#> [1] "2025-11-28 09:15:41 GMT" "2025-11-28 00:00:00 GMT"
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

## Casting and Coercion

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
#> 1 gcd(x)        900ns    1.3µs   672595.        0B     67.3
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       29.8ms   35.5ms      27.1        0B        0
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
#> 1 x * 10 * 20 + 1 - 1                 1.42ms 2.02ms      447.    7.63MB     36.4
#> 2 set_subtract(set_add(set_multiply(… 3.14ms 3.66ms      247.        0B      0
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
#> 1 do.call(c, x)   2.75ms   3.45ms      269.     781KB   261.  
#> 2 c_(.args = x)   1.03ms   1.14ms      796.     781KB     4.16

# Matches the speed of `unlist()` without removing attributes
unlist(list(Sys.Date()), recursive = FALSE)
#> [1] 20420
c_(.args = list(Sys.Date()))
#> [1] "2025-11-28"
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
#> 1 shallow_copy(iris)    400ns    500ns  1596271.    6.34KB        0
mark(deep_copy(iris))
#> # A tibble: 1 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 deep_copy(iris)    700ns    1.6µs   441902.    9.34KB     44.2
mark(semi_copy(iris))
#> # A tibble: 1 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 semi_copy(iris)    700ns      2µs   342563.    9.36KB        0
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
#> 1 semi_copy(df)    117µs    140µs     5437.    3.81MB     234.
#> 2 deep_copy(df)    235µs    328µs     2056.    7.63MB     207.
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
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 add_length_class(integer(10^6))      874µs  874µs     1144.    3.81MB        0
#> 2 add_length_class_in_place(integer(1… 339µs  339µs     2952.    3.81MB        0
mark(
  add_length_class(integer(10^6)),
  add_length_class_in_place(integer(10^6)),
  iterations = 1
)
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 add_length_class(integer(10^6))      592µs  592µs     1690.    3.81MB        0
#> 2 add_length_class_in_place(integer(1… 396µs  396µs     2525.    3.81MB        0
  

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
#> 1 cheapr_which   2.49ms   2.84ms      287.    3.82MB     15.7
#> 2 base_which    564.8µs  694.2µs     1138.    7.63MB    106.
x <- rep(FALSE, 10^6)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    118µs    125µs     6766.        0B       0 
#> 2 base_which      225µs    242µs     3435.    3.81MB     107.
x <- c(rep(TRUE, 5e05), rep(FALSE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   1.24ms   1.61ms      518.    1.91MB     8.60
#> 2 base_which    521.1µs  790.7µs     1144.    7.63MB    77.5
x <- c(rep(FALSE, 5e05), rep(TRUE, 1e06))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which    893µs   1.26ms      735.    3.81MB     26.1
#> 2 base_which      721µs   1.01ms      866.    9.54MB    101.
x <- sample(c(TRUE, FALSE), 10^6, TRUE)
x[sample.int(10^6, 10^4)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which  596.1µs  679.9µs     1319.    1.89MB     27.1
#> 2 base_which     3.72ms   4.19ms      220.     5.7MB     13.1
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
#> 1 cheapr_factor   9.05ms   10.3ms     96.9     6.13MB     4.40
#> 2 base_factor   318.69ms  318.7ms      3.14   27.84MB     3.14
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   3.12ms   4.29ms    221.      1.55MB     4.41
#> 2 base_factor   625.42ms 625.42ms      1.60   22.79MB     0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor  78.06ms  96.82ms    10.3      17.4MB    1.72 
#> 2 base_factor      3.42s    3.42s     0.292    54.4MB    0.292
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   5.03ms   9.29ms     110.     3.49MB     7.31
#> 2 base_factor    49.14ms  52.96ms      18.9   39.89MB    47.2
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
#> 1 cheapr_intersect    2.4ms   3.81ms      257.    1.55MB     4.55
#> 2 base_intersect     4.67ms   7.58ms      130.    6.41MB     7.49
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff   2.11ms   3.06ms      313.    2.15MB     7.00
#> 2 base_setdiff     4.34ms   5.12ms      186.    6.96MB    12.7
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.31ms    1.5ms      605.  781.34KB     6.77
#> 2 base         2.22ms   2.58ms      360.    2.53MB     9.53
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       1.22ms   1.36ms      670.  792.29KB     6.89
#> 2 base         2.28ms   2.67ms      361.    2.91MB     9.51
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
#> 1 cheapr_cut   14.9ms   17.3ms      57.1    3.92MB     4.57
#> 2 base_cut     25.6ms   32.9ms      30.9   15.32MB     7.72
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
#> 1 if_else_(x >= 0, 1, -1)            2.77ms  3.38ms     286.     11.4MB     56.3
#> 2 ifelse(x >= 0, 1, -1)             17.86ms 20.22ms      49.5    53.4MB    158. 
#> 3 data.table::fifelse(x >= 0, 1, -…  5.81ms  6.91ms     141.     11.4MB     17.3
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
#> 1 "case(x >= 0 ~ \"pos\", x < 0 ~ \"…   19ms   22ms      45.7    28.8MB     22.8
#> 2 "data.table::fcase(x >= 0, \"pos\"… 16.2ms 18.1ms      54.4    26.7MB     32.7
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
#> 1 val_match(x, 1 ~ Inf, 2 ~ -Inf, …  3.59ms  4.79ms     205.     8.79MB     27.5
#> 2 case(x == 1 ~ Inf, x == 2 ~ -Inf… 14.42ms 18.71ms      52.8   27.63MB     35.2
#> 3 data.table::fcase(x == 1, Inf, x… 11.12ms 13.67ms      66.8   30.52MB     61.7
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
#> 1 get_breaks(x, 20)   61.5µs   67.8µs    12348.        0B      0  
#> 2 pretty(x, 20)      399.5µs  614.8µs     1541.    1.91MB     55.9

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

<img src="man/figures/README-unnamed-chunk-60-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-60-2.png" width="100%" />
