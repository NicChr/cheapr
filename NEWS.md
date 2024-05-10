# cheapr 0.9.2

* A signed integer overflow bug in `lag2_` has been fixed. This occurred when 
supplying `NA` lags. 

* `lag2_` no longer fills the names of named vectors when the `fill` value is supplied.

# cheapr 0.9.1 (05-May-2024)

* New function `recycle` to help recycle R objects to a common size.

* The `set` functions that update by reference are now ALTREP aware and
take a copy when the input is an ALTREP object.

* New function `lag2_` as a generalised solution for complex lags. It supports
dynamic lag vectors, lags using an order vector, and custom run lengths. 
It doesn't support updating by reference or long vectors.

# cheapr 0.9.0 (22-Apr-2024)

* New function `lag_` for very fast lags and leads on vectors and data frames.
It includes a `set` argument allowing users to create a lagged vector 
by reference without copies.

* `set_round` has been amended to improve floating point accuracy.

# cheapr 0.8.0 (12-Apr-2024)

* New 'set' Math operations inspired by 'data.table' and 'collapse' 
that transform data by reference.

* Fixed an inconsistency of when `sequence_()` would error when supplied with 
a zero-length size argument.

* Fixed a protection stack imbalance in `count_val(x)` when `x` is `NULL`.

* `sset` has been optimised for wide data frames with many variables. 
It is also faster when applied to a data frame with dates, date-times and factors.

* In `sset`, when `i` is a logical vector it must match the length of x.

* `sset` can now handle 'ALTREP' compact real sequences as well.

# cheapr 0.5.0 (5-Apr-2024)

* `sset` is now parallelised when `i` is an 'ALTREP'
compact integer sequence, e.g. `sset(x, 1:10)`.

* `sset` now has an internal range-based subset method for 
'ALTREP' integer sequences made using  `:` for example.

* New function `count_val` as a cheaper alternative to e.g. `sum(x == val)`.

* Negative indexing in `sset` has been improved. 
It is also now partially parallelised.

* Setting `recursive` to false should now be faster.

* 'overview' objects gain an additional list element "print_digits" which 
is passed to the print method in order to correctly round the summary statistics 
without affecting the 'cheapr.digits' option globally.

* `factor_` and `na_rm` now handle data frames.

* A bug in `sset.data.table` that caused further set calculations to produce 
warnings has been fixed.

* `is_na.POSIXlt` and `sset.POSIXlt` have been rewritten to handle unbalanced 
'POSIXlt' objects.

# cheapr 0.4.0 (25-Mar-2024)

* New function `sset` to consistently subset data frame rows and vectors in 
general.

* `overview` now always returns an object of class "overview". It also returns
the number of observations instead of rows so that it makes sense 
for vector summaries as well as data frame summaries.

* `sequence_` has been optimised and rewritten in C++. It now only checks for
integer overflow when both `from` and `by` are integer vectors.

* The internal function `list_as_df` has been rewritten in C++.

# cheapr 0.3.0 (18-Mar-2024)

* New function `overview` as a cheaper alternative to `summary`.

* All of the `NA` handling functions now fall back to using `is.na` if an appropriate
method cannot be found.

* More support has been added for all objects with an `is.na` method.

# cheapr 0.2.0 (06-Mar-2024)

* `is_na` has been added as an S3 generic function which is parallelised and  internally falls back
on `is.na` if there are no suitable methods.

* Additional list utility functions have been added.

* Limited support for `vctrs_rcrd` objects has been added again. 

* `num_na` and similar functions no longer treat empty data frame rows as single observations but instead return the total number of `NA` values in the data frame.

* Fixed a bug in `row_na_counts` and `col_na_counts` that would cause the 
session to crash when a column variable was a list.

* For the time being, vctrs 'vctrs_rcrd' objects are no longer supported though
this support may be re-added in the future.

# cheapr 0.1.0 (05-Mar-2024)

* CRAN submission accepted.
