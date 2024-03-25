# cheapr 0.4.0

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
