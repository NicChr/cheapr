#' Count the number of occurrences of a value.
#'
#' @description
#' This is a programmer's version of `sum(x == value)` to count the number of
#' occurrences of a value without creating a potentially large logical vector.
#'
#' @param x A vector, list, data frame or matrix.
#' @param value A value with which to count the frequency of.
#' @param recursive Should the function be applied recursively to lists?
#'
#' @details
#' This is a generalisation of `num_na()` and as such the identity
#' `count_val(x, NA) == num_na(x)` will always hold.
#'
#' @returns
#' A count of the number of times `value` appears in `x`.
#'
#' @export
count_val <- function(x, value, recursive = TRUE){
  cpp_count_val(x, value, recursive)
}
