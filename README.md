
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

## cheapr C/C++ API

cheapr includes a feature-rich C API meant as an alternative to the R C
API. Features include efficient vector creation, replacement, subsetting
and comparison. The cheapr C API also includes tools to efficiently work
with C types while carefully handling `NA` values and other R-specific
values.  
  
**Important note:** Despite the API requiring C++17 and making heavy use
of C++ templates, C++ objects cannot be safely used. More specifically,
C++ objects with destructors cannot be safely used because when R
encounters an error, it terminates in a way that doesn’t run these
destructors, leading to memory leaks. This applies to any C++ object
with a destructor like `std::string` and `std::vector`.

We plan on creating a C++ safe API in the future, but if you still
really want to use C++ objects (like `std::vector`) along with the
cheapr or R C API then you must call any C API functions with either
`cpp11::safe[]` or cheapr’s `r_safe()`. For more info see cpp11’s FAQ
section 15 [Should I call cpp11::unwind_protect()
manually?](https://cpp11.r-lib.org/articles/FAQ.html)

The reason we have opted for a mostly C API with C++ templates is for
maximal speed, efficiency and interoperability with the R C API.

We aim to create a fully C++ compatible API in an upcoming release.

## The double-edged `SEXP`

R internally labels nearly all allocated objects as `SEXP`
(S-Expression). While this provides much flexibility in inputs and
outputs, it also creates ambiguity regarding the underlying type: a
`SEXP` might represent a vector (like `INTSXP`), a single element (like
a string `CHARSXP`), or a language object (like `SYMSXP`). In the
standard R C API, these can be disambiguated at runtime using
`TYPEOF()`, which makes meta-programming cumbersome and difficult.  
Furthermore because of the lack of generic functions in the R C API,
dedicated ones must be used that are specific to the type of object one
is working with.

`cheapr` solves this by introducing distinct C++ types for R’s
non-vector elements, enabling optimised template specialisations. The
underlying types are C-based and do not need C++ destructors.

## `cheapr` C/C++ Types

The following types should be preferred over standard R API types to
ensure cheapr templates function correctly.

- **`r_bool_t`** – Represents R’s tri-state logical (TRUE, FALSE, NA).
  Is internally `int`

- **`r_string_t`** – Alias for `CHARSXP` (internal R strings)

- **`r_symbol_t`** – Alias for `SYMSXP` (R symbols).

- **`r_byte_t`** – Alias for `Rbyte` (raw elements).

- **`r_complex_t`** – An alternative to `Rcomplex`

## Disambiguating `CHARSXP` and `SYMSXP`

By introducing `r_string_t` and `r_symbol_t`, we remove the ambiguity of
`SEXP` in template meta-programming in the context of vectors.

In `cheapr`, the `SEXP` type is effectively reserved for **list
elements** (or generic R objects).

- `new_vector<int>` creates an integer vector (`INTSXP`)

- `new_vector<r_string_t>` creates a character vector (`STRSXP`)

- `new_vector<SEXP>` creates a **list** (`VECSXP`)

Since a list is the only vector capable of holding other R objects,
`new_vector<SEXP>` unambiguously resolves to creating a new list.

## `r_cast<>`

Standard C++ casts (`static_cast<>`) are unaware of R’s special `NA`
values.

- `static_cast<int>(NA_REAL)` results in `INT_MIN` not `NA_INTEGER`.

`cheapr` provides **`r_cast<>`**, a powerful cast that handles R
specific values correctly:

- `r_cast<int>(na::real)` returns `na::integer`.

- `r_cast<double>(na::integer)` returns `na::real`.

- `r_cast<int>(limits::r_pos_inf)` returns `na::integer`

- `r_cast<r_symbol_t>("my_string")` returns the R symbol `my_string`

- `r_cast<r_string_t>(R_ClassSymbol)` returns the internal R string
  “class”

## `NA` constants

All R-specific `NA` values are available in the `na` namespace

| `cheapr` Constant | C/C++ Type    | R C API Macro  |
|:------------------|:--------------|:---------------|
| `na::logical`     | `r_bool_t`    | `NA_LOGICAL`   |
| `na::integer`     | `int`         | `NA_INTEGER`   |
| `na::real`        | `double`      | `NA_REAL`      |
| `na::complex`     | `r_complex_t` | Not applicable |
| `na::string`      | `r_string_t`  | `NA_STRING`    |
| `na::nil`         | `SEXP`        | `R_NilValue`   |

## Creating Vectors

Use `new_vector<T>(size, fill_value)` to create R vectors efficiently.

``` cpp
// Create integer vector of length 10
SEXP v1 = new_vector<int>(10);

// Create character vector of length 20
SEXP v2 = new_vector<r_string_t>(20);

// Create integer vector filled with NAs
SEXP v3 = new_vector<int>(100, na::integer);
```

**Note:** Large vectors (except character/list vectors) are initialised
using multiple threads which be adjusted via `cheapr::set_threads()`

## Creating Scalars

To create length-one vectors (scalars) efficiently:

``` cpp
// Create scalar character vector
SEXP s1 = as_vector("hello");
// Similarly
SEXP s2 = as_vector(r_cast<r_string_t>("hello"));
```

Alternatively, `r_cast<SEXP>` can convert any object to a SEXP. Note
that `r_string_t` are always converted to length-one character vectors
and `r_symbol_t` objects are always converted to length-one lists.

## Combining Values

`cheapr::combine()` (the backend for `c_()`) unifies C++ types and R
vectors into a single vector. It automatically determines the common
type (see `?cast`).

``` cpp
// Returns a character vector: c("1", "2.5", "str1")
combine(1, 2.5, "str1");
```

## Creating Lists

Lists are created via `make_list()` with optional named arguments:

``` cpp
// Equivalent to list(1, 2, 3)
SEXP l1 = make_list(1, 2, 3);

// Equivalent to list(x = 1, y = 2, z = 3)
SEXP l2 = make_list(arg("x") = 1, arg("y") = 2, arg("z") = 3);
```

## Math Utilities

The `cheapr::math` namespace provides safe, NA aware implementations of
common math functions, including `math::r_abs()`, `math::r_round()`
(round-to-even), `math::r_gcd()` and more.

### Getting & setting vector values

In the R C API you would typically either use a pointer to the
underlying data (very efficient) to get/set data, or you would use a
function like `INTEGER_ELT` or `SET_INTEGER_ELT` to respectively get and
set an int value (less efficient).

cheapr uses overloaded template functions `get_value<>` and
`set_value<>` to simplify this without sacrificing efficiency.

``` cpp
SEXP int_vec = SHIELD(as_vector(1)); // SHIELD is the same as Rf_protect

// Get the vector's first value

// Have to specify the type int in get_value<>
int first_val = get_value<int>(int_vec, 0)
Rprintf("first value: %d", first_val);

// Set the vector's first value

// `set_value<>` can deduce that it is an int based on the 
// type of value being set
set_value(int_vec, 0, 42);
// or more explicitly 
set_value<int>(int_vec, 0, 42);
Rprintf("first value: %d", get_value<int>(int_vec, 0));
```

``` r
first value: 1
first value: 42
```

The getters/setters can also accept pointers directly

``` cpp
auto *vec_ptr = vector_ptr<int>(int_vec);

set_value(vec_ptr, 0, 123);
auto first_val = get_value<int>(vec_ptr, 0);

Rprintf("\nfirst value: %d", first_val);
```

``` r
first value: 123
```

Don’t forget to unprotect all previously protected objects

``` cpp
YIELD(1); // Same as Rf_unprotect
```

Some notes:

- Strings are converted to UT8 almost everywhere

- Coercion between the R API and C API may not always match, especially
  for data frames. This is something that is likely to be improved in
  the future.

- Because cheapr tries to treat data frames as vectors wherever
  possible, functions like `vec::length()` return the number of rows of
  a data frame. Use `vec::old_length()` to match `Rf_xlength()`

### An alternative to switch(TYPEOF())

Sometimes you need to apply different code to handle the different
nuances of your vector type. Let’s showcase an example of this by
calculating the sum of a vector using plain R C API style code.

``` cpp
[[cpp11::register]]
double r_style_sum(SEXP x){
  double sum = 0;

  R_xlen_t n = vec::old_length(x);

  switch (TYPEOF(x)){
  case LGLSXP: {
    int *p_x = LOGICAL(x);
    for (R_xlen_t i = 0; i < n; ++i){
      sum += p_x[i] == NA_LOGICAL ? NA_REAL : static_cast<double>(p_x[i]);
    }
    break;
  }
  case INTSXP: {
    int *p_x = INTEGER(x);
    for (R_xlen_t i = 0; i < n; ++i){
      sum += p_x[i] == NA_INTEGER ? NA_REAL : static_cast<double>(p_x[i]);
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    for (R_xlen_t i = 0; i < n; ++i){
      sum += p_x[i];
    }
    break;
  }
  default: {
    Rf_error("`x` must be an integer or double vector in `r_style_sum`");
  }
  }
  return sum;
}
```

Notice that every branch has slightly different syntax. You could reduce
the first two branches into one by using `INTEGER()` for logical vector
inputs but I digress.

In cheapr this can be done in a more generic way using a C++ lambda and
cheapr’s templated functions like `r_cast<>`

``` cpp
[[cpp11::register]]
double cheapr_sum(SEXP x){

  double out = visit_vector(x, [&](auto p_x) -> double {

    // Data type of vector elements
    using data_t = std::decay_t<std::remove_pointer_t<decltype(p_x)>>;

    double sum = 0;

    // This includes C++ arithmetic types and `r_bool_t`
    // constexpr allows this if statement to be evaluated at compile-time
    if constexpr (is_r_arithmetic_v<data_t>){

      R_xlen_t n = vec::old_length(x);

      for (R_xlen_t i = 0; i < n; ++i){
        // r_cast will naturally handle NA values
        sum += r_cast<double>(p_x[i]);
      }
    } else {
      Rf_error("`x` must be numeric in `cheapr_sum`");
    }
    return sum;
  });

  return out;
}
```

## cheapr R API

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

Let’s first load the required packages

``` r
library(cheapr) 
library(bench)
```

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
#> 1 row_na_counts(m)   483.9µs  583.6µs     1618.   13.08KB      0  
#> 2 rowSums(is.na(m))   2.99ms   3.57ms      281.    3.84MB     25.9
# Number of NA values by col
mark(col_na_counts(m), 
     colSums(is.na(m)))
#> # A tibble: 2 × 6
#>   expression             min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 col_na_counts(m)    1.36ms   1.57ms      603.   13.08KB      0  
#> 2 colSums(is.na(m))   1.27ms   1.54ms      598.    3.82MB     60.0
```

`is_na` is a multi-threaded alternative to `is.na`

``` r
x <- rnorm(10^6) |> 
  na_insert(10^5)
set_threads(4)
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      596µs    722µs     1244.    3.81MB     228.
#> 2 is_na(x)      161µs    250µs     4154.    3.82MB     414.
set_threads(1)
mark(is.na(x), is_na(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 is.na(x)      592µs    730µs     1277.    3.81MB     134.
#> 2 is_na(x)      390µs    479µs     2001.    3.81MB     200.
set_threads(2)

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
#> 1 is.na(hours)    1.18s    1.18s     0.849   61.05MB     1.70
#> 2 is_na(hours)   3.79ms   4.69ms   194.       7.65MB    13.9
```

It differs in 2 regards:

- List vectors are never regarded as `NA` even if they contain `NA`
  values.
- For data frames, `is_na` returns a logical vector where `TRUE` defines
  an empty row of only `NA` values.

``` r
# List example
is.na(list(NA, NULL, 10))
#> [1]  TRUE FALSE FALSE
is_na(list(NA, NULL, 10))
#> [1] FALSE FALSE FALSE

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
#> 1 overview(df, hist = FALSE)   72.4ms   79.7ms      12.6      512B        0
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
#> 1 sset(x, x %in_% y)   88.1µs    109µs     7769.      86KB     8.32
#> 2 sset(x, x %in% y)   132.3µs    160µs     5503.     286KB    19.2 
#> 3 x[x %in% y]         137.8µs    208µs     4424.     325KB    15.1
```

`sset` uses an internal range-based subset when `i` is an ALTREP integer
sequence of the form m:n.

``` r
mark(sset(df, 0:10^5), df[0:10^5, , drop = FALSE])
#> # A tibble: 2 × 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, 0:10^5)            462.6µs  607.1µs     1484.    3.43MB     68.2
#> 2 df[0:10^5, , drop = FALSE]   6.26ms    7.2ms      128.    4.83MB     13.5
```

It also accepts negative indexes

``` r
mark(sset(df, -10^4:0), 
     df[-10^4:0, , drop = FALSE],
     check = FALSE) # The only difference is the row names
#> # A tibble: 2 × 6
#>   expression                       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sset(df, -10^4:0)             2.66ms   4.38ms     233.     18.9MB     110.
#> 2 df[-10^4:0, , drop = FALSE]  24.48ms  24.48ms      40.8    72.5MB     858.
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
#> [1] "2025-12-25" "2025-12-25"
#> [1] "2025-12-25 12:08:39 GMT" "2025-12-25 00:00:00 GMT"
c_(today, now);c_(now, today) # cheapr
#> [1] "2025-12-25 00:00:00 GMT" "2025-12-25 12:08:39 GMT"
#> [1] "2025-12-25 12:08:39 GMT" "2025-12-25 00:00:00 GMT"
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
data.table::setDTthreads(2);
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
#> 1 gcd(x)        800ns      1µs   791755.        0B     79.2
x <- seq(0, 10^6, 0.5)
mark(gcd(x))
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 gcd(x)       31.3ms   34.8ms      27.5        0B        0
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

## Parallel math

Many common math functions can be found in cheapr and are parallelised

``` r
set_threads(4)

x <- abs(rnorm(3e05))

mark(negate_(x), -x)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 negate_(x)   87.5µs    154µs     5236.    2.29MB    112. 
#> 2 -x            260µs    396µs     2274.    2.29MB     41.2
mark(abs_(x), abs(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 abs_(x)      98.1µs    146µs     5549.    2.29MB    109. 
#> 2 abs(x)      179.9µs    233µs     3011.    2.29MB     60.6
mark(floor_(x), floor(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 floor_(x)   114.3µs  189.8µs     4526.    2.29MB     93.2
#> 2 floor(x)     1.48ms   1.71ms      588.    2.29MB     10.9
mark(ceiling_(x), ceiling(x))
#> # A tibble: 2 × 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 ceiling_(x)  113.4µs    192µs     4519.    2.29MB     90.7
#> 2 ceiling(x)    1.29ms   1.55ms      593.    2.29MB     13.3
mark(trunc_(x), trunc(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 trunc_(x)   110.5µs 177.15µs     5031.    2.29MB    108. 
#> 2 trunc(x)      1.3ms   1.49ms      575.    2.29MB     10.9
mark(sqrt_(x), sqrt(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 sqrt_(x)      137µs    245µs     4040.    2.29MB     90.2
#> 2 sqrt(x)       636µs    671µs     1261.    2.29MB     27.9
mark(log_(x), log(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 log_(x)      1.76ms   2.16ms      437.     2.3MB     9.55
#> 2 log(x)       6.91ms   7.16ms      123.    2.29MB     2.05
mark(exp_(x), exp(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 exp_(x)      3.42ms   3.96ms     227.     2.29MB     6.30
#> 2 exp(x)      11.82ms  12.38ms      72.1    2.29MB     0
mark(round_(x), round(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 round_(x)    1.08ms   1.33ms      683.     2.3MB    15.4 
#> 2 round(x)     8.69ms   9.07ms      100.    2.29MB     2.00
mark(signif_(x), signif(x))
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 signif_(x)    3.1ms      4ms     240.      2.3MB     6.36
#> 2 signif(x)    14.5ms   14.8ms      60.4    2.29MB     0
mark(add_(x, x), x + x)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 add_(x, x)   97.7µs    152µs     4507.    2.29MB     98.6
#> 2 x + x       184.3µs    215µs     3250.    2.29MB     69.5
mark(subtract_(x, x), x - x)
#> # A tibble: 2 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 subtract_(x, x)   98.1µs    149µs     5439.    2.29MB    118. 
#> 2 x - x            170.5µs    254µs     3137.    2.29MB     66.2
mark(multiply_(x, x), x * x)
#> # A tibble: 2 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 multiply_(x, x)   92.2µs    130µs     5589.    2.29MB    114. 
#> 2 x * x            179.3µs    320µs     2647.    2.29MB     54.3
mark(divide_(x, x), x / x)
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 divide_(x, x)   90.2µs    228µs     3963.    2.29MB     86.1
#> 2 x/x            191.1µs    225µs     2684.    2.29MB     57.1
mark(pow_(x, 2), x^2)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 pow_(x, 2)    120µs    159µs     4308.    2.29MB     87.5
#> 2 x^2           181µs    224µs     3064.    2.29MB     57.5

set_threads(2)
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
#> 1 x * 10 * 20 + 1 - 1                 1.39ms 2.25ms      487.    7.63MB     41.0
#> 2 set_subtract(set_add(set_multiply(… 1.87ms 1.97ms      473.        0B      0
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
#> 1 do.call(c, x)   2.72ms   3.91ms      238.     781KB   133.  
#> 2 c_(.args = x)   1.02ms   1.15ms      773.     781KB     4.15

# Matches the speed of `unlist()` without removing attributes
unlist(list(Sys.Date()), recursive = FALSE)
#> [1] 20447
c_(.args = list(Sys.Date()))
#> [1] "2025-12-25"
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
#> 1 shallow_copy(iris)    300ns    400ns  2070522.    6.34KB        0
mark(deep_copy(iris))
#> # A tibble: 1 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 deep_copy(iris)    3.6µs    5.6µs   152729.    10.6KB     15.3
mark(semi_copy(iris))
#> # A tibble: 1 × 6
#>   expression           min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 semi_copy(iris)    2.5µs    4.5µs   161579.    9.36KB     16.2
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
#> 1 semi_copy(df)  211.4µs  249.2µs     3267.    3.81MB     160.
#> 2 deep_copy(df)    1.4ms   1.69ms      500.   15.26MB     125.
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
  add_length_class(integer(10^5)),
  add_length_class_in_place(integer(10^5)),
  iterations = 1
)
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 add_length_class(integer(10^5))      175µs  175µs     5711.     391KB        0
#> 2 add_length_class_in_place(integer(1… 166µs  166µs     6013.     391KB        0
mark(
  add_length_class(integer(10^5)),
  add_length_class_in_place(integer(10^5)),
  iterations = 1
)
#> # A tibble: 2 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 add_length_class(integer(10^5))      117µs  117µs     8562.     391KB        0
#> 2 add_length_class_in_place(integer(1… 144µs  144µs     6920.     391KB        0
  

# R detected that the vector we created had been modified (because it was)
# and created a copy
# When we add the attributes in-place to our fresh object, no copies are
# made
```

## ‘Cheaper’ Base R alternatives

### which

``` r
x <- rep(TRUE, 10^5)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   82.2µs    240µs     2683.     397KB     10.5
#> 2 base_which     69.8µs   79.9µs     8347.     781KB     56.9
x <- rep(FALSE, 10^5)
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   18.9µs   33.3µs    27903.        0B      0  
#> 2 base_which     23.3µs   26.4µs    29181.     391KB     93.7
x <- c(rep(TRUE, 5e04), rep(FALSE, 1e05))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   53.2µs   89.3µs    10800.     195KB     15.6
#> 2 base_which     57.7µs    124µs     9080.     781KB     54.6
x <- c(rep(FALSE, 5e04), rep(TRUE, 1e05))
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which  100.5µs    165µs     5699.     391KB     18.6
#> 2 base_which     77.2µs    132µs     6959.     977KB     53.9
x <- sample(c(TRUE, FALSE), 10^5, TRUE)
x[sample.int(10^5, 10^3)] <- NA
mark(cheapr_which = which_(x),
     base_which = which(x))
#> # A tibble: 2 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_which   69.1µs    108µs     8900.     193KB    13.2 
#> 2 base_which    367.5µs    410µs     2080.     584KB     8.43
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
#> 1 cheapr_factor   9.19ms   10.7ms     78.4     6.13MB     6.28
#> 2 base_factor   313.54ms  350.1ms      2.86   27.84MB     0
mark(cheapr_factor = factor_(x, order = FALSE), 
     base_factor = factor(x, levels = unique(x)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   2.96ms   3.28ms    253.      1.55MB     4.29
#> 2 base_factor    615.9ms  615.9ms      1.62   22.79MB     0
mark(cheapr_factor = factor_(y), 
     base_factor = factor(y))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor  66.29ms  85.36ms    10.1      17.4MB    3.38 
#> 2 base_factor      3.63s    3.63s     0.275    54.4MB    0.275
mark(cheapr_factor = factor_(y, order = FALSE), 
     base_factor = factor(y, levels = unique(y)))
#> # A tibble: 2 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_factor   4.38ms   5.55ms     121.     3.49MB     5.37
#> 2 base_factor    46.78ms   54.4ms      17.5   39.89MB    14.0
```

### intersect & setdiff

``` r
x <- sample.int(10^6, 10^4, TRUE)
y <- sample.int(10^6, 10^4, TRUE)
mark(cheapr_intersect = intersect_(x, y, dups = FALSE),
     base_intersect = intersect(x, y))
#> # A tibble: 2 × 6
#>   expression            min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>       <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_intersect    139µs    217µs     4066.     167KB     9.67
#> 2 base_intersect      298µs    343µs     2387.     736KB    28.4
mark(cheapr_setdiff = setdiff_(x, y, dups = FALSE),
     base_setdiff = setdiff(x, y))
#> # A tibble: 2 × 6
#>   expression          min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>     <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr_setdiff    147µs    219µs     3972.     254KB     9.97
#> 2 base_setdiff      313µs    373µs     2197.     812KB    18.1
```

### `%in_%` and `%!in_%`

``` r
mark(cheapr = x %in_% y,
     base = x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       79.2µs    102µs     8163.    78.2KB     6.67
#> 2 base          122µs    151µs     5209.   284.5KB    14.6
mark(cheapr = x %!in_% y,
     base = !x %in% y)
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 cheapr       74.3µs   97.7µs     9077.    89.2KB     6.71
#> 2 base        133.5µs  157.7µs     4874.   323.6KB    14.7
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
#> 1 cheapr_cut   13.6ms     15ms      62.4    3.92MB     4.62
#> 2 base_cut     25.4ms   28.3ms      32.0   15.32MB     7.39
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
#> 1 if_else_(x >= 0, 1, -1)            2.58ms  3.09ms     276.     11.4MB     65.8
#> 2 ifelse(x >= 0, 1, -1)             17.75ms 19.03ms      52.9    53.4MB    172. 
#> 3 data.table::fifelse(x >= 0, 1, -…  3.91ms  4.92ms     197.     11.4MB     35.6
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
#> 1 "case(x >= 0 ~ \"pos\", x < 0 ~ \"… 17.6ms 19.3ms      48.3    28.8MB     48.3
#> 2 "data.table::fcase(x >= 0, \"pos\"… 15.6ms 18.7ms      50.7    26.7MB     32.6
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
#> 1 val_match(x, 1 ~ Inf, 2 ~ -Inf, …  3.14ms  3.57ms     250.     8.78MB     48.3
#> 2 case(x == 1 ~ Inf, x == 2 ~ -Inf… 13.46ms 17.41ms      57.2   27.62MB     52.0
#> 3 data.table::fcase(x == 1, Inf, x… 10.97ms 12.78ms      71.2   30.52MB     35.6
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
#> 1 get_breaks(x, 20)   61.5µs   65.6µs    11794.        0B      0  
#> 2 pretty(x, 20)      414.2µs  517.5µs     1671.    1.91MB     24.6

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

<img src="man/figures/README-unnamed-chunk-53-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-53-2.png" width="100%" />
