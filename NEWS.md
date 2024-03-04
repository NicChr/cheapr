# cheapr (Development version)

* `num_na` and similar functions no longer treat empty data frame rows as single observations but instead return the total number of `NA` values in the data frame.

* Fixed a bug in `row_na_counts` and `col_na_counts` that would cause the 
session to crash when a column variable was a list.

* For the time being, vctrs 'vctrs_rcrd' objects are no longer supported though
this support may be re-added in the future.

# cheapr 0.1.0 (Pending)

* CRAN submission.
