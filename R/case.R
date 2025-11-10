#' A cheapr case-when and switch
#'
#'
#' @description
#' `case` and `val_match` are cheaper alternatives to `dplyr::case_when` and
#' `dplyr::case_match` respectively.
#'
#'
#' @param ... Logical expressions or scalar values in the case of `val_match`.
#' @param .x Vector used to switch values.
#' @param .default Catch-all value or vector.
#'
#' @seealso [if_else_]
#'
#' @returns
#' A vector the same length as `.x` or same length as the first condition in the
#' case of `case`, unless the condition length is smaller than the
#' rhs, in which case the length of the rhs is used.
#'
#' @details
#' `val_match()` is a very efficient special case of the
#' `case()` function when all lhs expressions are scalars,
#' i.e. length-1 vectors. RHS expressions can be vectors the
#' same length as `.x`.
#' The below 2 expressions are equivalent.
#'
#' \preformatted{
#' val_match(
#'   x,
#'   1 ~ "one",
#'   2 ~ "two",
#'   .default = "Unknown"
#'  )
#' case(
#'   x == 1 ~ "one",
#'   x == 2 ~ "two",
#'   .default = "Unknown"
#'  )
#' }
#'
#'
#' @rdname case
#' @export
case <- function(..., .default = NULL){
  conditions <- exprs(...)
  force(.default)
  n_conditions <- length(conditions)

  if (n_conditions == 0){
    stop("Please supply >= 1 case when statements")
  }

  rhs_list <- new_list(n_conditions)
  lgl_list <- new_list(n_conditions)

  # Evaluate first condition

  condition <- conditions[[1L]]
  if (as.character(condition[[1L]]) != "~" || length(condition) != 3){
    stop("Please supply formula expressions")
  }

  # Deal with RHS values
  for (i in seq_along(conditions)){
    rhs_list[[i]] <- eval(conditions[[i]][[3]], envir = parent.frame())
  }

  if (is.null(.default)){
    rhs_list <- cast_common(.args = rhs_list)
  } else {
    rhs_common <- cast_common(.args = c_(rhs_list, list(.default)))
    rhs_list <- rhs_common[seq_along(rhs_list)]
    .default <- rhs_common[[length(rhs_common)]]
  }

  lgl <- eval(condition[[2L]], envir = parent.frame())
  if (!is.logical(lgl)){
    stop("All LHS expressions must be logical vectors")
  }

  rhs <- rhs_list[[1]]
  if (length(lgl) == 0){
    return(sset(rhs, 0))
  }

  if (length(lgl) < vector_length(rhs)){
    lgl <- rep_len_(lgl, vector_length(rhs))
  }
  if (!vector_length(rhs) %in% c(1, length(lgl))){
    stop("rhs must be of length 1 or length of lhs")
  }

  out_size <- length(lgl)

  true_locs <- val_find(lgl, TRUE)

  out <- na_init(rhs, out_size)

  if (length(rhs) == 1){
    out <- .Call(`_cheapr_cpp_replace`, out, true_locs, rhs, FALSE, FALSE)
  } else {
    out <- .Call(`_cheapr_cpp_replace`, out, true_locs, sset(rhs, true_locs), FALSE, FALSE)
  }

  # Initialise a logical that will keep track of
  # which conditions have been met
  lgl_or <- deep_copy(lgl)

  if (n_conditions >= 2){
    lgl3 <- logical(length(lgl))
    for (i in 2:n_conditions){

      condition <- conditions[[i]]
      if (as.character(condition[[1L]]) != "~" || length(condition) != 3){
        stop("Please supply formula expressions")
      }
      rhs <- rhs_list[[i]]
      lgl <- eval(condition[[2L]], envir = parent.frame())

      if (!is.logical(lgl)){
        stop("All LHS expressions must be logical vectors")
      }

      lgl <- rep_len_(lgl, out_size)
      if (!vector_length(rhs) %in% c(1, out_size)){
        stop("rhs must be of length 1 or length of lhs")
      }

      # 1. Copy elements from user condition into our temporary lgl3 vector
      # 2. Replace previous evaluated `TRUE` elements with `FALSE`
      # to keep priority of former conditions
      # 3. Do `lgl_or | lgl` to keep track of evaluated `TRUE` elements

      cpp_set_copy_elements(source = lgl, target = lgl3)
      replace_(lgl3, val_find(lgl_or, TRUE), FALSE, in_place = TRUE, quiet = TRUE)
      cpp_set_or(lgl_or, lgl)

      true_locs <- val_find(lgl3, TRUE)
      if (length(true_locs) > 0){
        if (vector_length(rhs) == 1){
          out <- .Call(`_cheapr_cpp_replace`, out, true_locs, rhs, FALSE, FALSE)
        } else {
          out <- .Call(`_cheapr_cpp_replace`, out, true_locs, sset(rhs, true_locs), FALSE, FALSE)
        }
      }
    }
  }

  if (!is.null(.default)){

    if (!vector_length(.default) %in% c(1, out_size)){
      stop("rhs must be of length 1 or length of lhs")
    }

    default_locs <- val_find(lgl_or, TRUE, invert = TRUE)

    if (vector_length(default_locs) > 0){
      if (vector_length(.default) == 1){
        out <- .Call(`_cheapr_cpp_replace`, out, default_locs, .default, FALSE, FALSE)
      } else {
        out <- .Call(`_cheapr_cpp_replace`, out, default_locs, sset(.default, default_locs), FALSE, FALSE)
      }
    }

  }
  out
}

#' @rdname case
#' @export
val_match <- function(.x, ..., .default = NULL){
  exprs <- exprs(...)

  N <- length(.x)
  n_exprs <- length(exprs)

  if (!(is.atomic(.x) || is.null(.x))){
    stop("`.x` must be an atomic vector")
  }

  if (!is.null(.default) && !length(.default) %in% c(1, N)){
    stop("`.default` must be length 1, `length(.x)` or `NULL`")
  }

  # Loop through formulas to check if all RHS expressions are
  # length-1 vectors, in which case we can use `fmatch()`
  # otherwise we use a slightly slower method

  # There shouldn't be too many expressions so looping twice
  # should carry minimal overhead

  rhs_all_scalars <- TRUE
  all_same_type <- TRUE
  # rhs_all_scalars <- is.null(.default) || length(.default) == 1

  lhs_list <- new_list(n_exprs)
  rhs_list <- new_list(n_exprs)

  for (i in seq_along(exprs)){
    expr <- exprs[[i]]
    if (as.character(expr[[1L]]) != "~" || length(expr) != 3){
      stop("Please supply formula expressions")
    }
    lhs <- eval(expr[[2L]], envir = parent.frame())
    if (length(lhs) != 1){
      stop("all lhs elements must be length 1 vectors")
    }
    rhs <- eval(expr[[3L]], envir = parent.frame())
    if (!length(rhs) %in% c(1, N)){
      stop("all rhs elements must be length 1 or `length(.x)`")
    }

    lhs_list[[i]] <- lhs
    rhs_list[[i]] <- rhs

    rhs_all_scalars <- rhs_all_scalars && length(rhs) == 1
    all_same_type <- all_same_type &&
      identical(typeof(.x), typeof(lhs))
  }

  # Are all LHS elements the same type as `x` and
  # are all RHS elements the same type?
  all_same_type <- all_same_type &&
    collapse::fnunique(
      vapply(if (is.null(.default)) rhs_list else c(rhs_list, list(.default)) , typeof, "",
             USE.NAMES = FALSE)
    ) == 1L

  if (rhs_all_scalars && (n_exprs > 3 || !all_same_type)){

    # Pre-allocate key-value pairs
      keys <- rep_len_(.x[NA_integer_], n_exprs)
      values <- keys

      # Assign keys and values
      for (i in seq_along(exprs)){
        keys[i] <- lhs_list[[i]]
        values[i] <- rhs_list[[i]]
      }

      # Add .default to values and ensure no-match locs point to that
      key_locs <- collapse::fmatch(.x, keys,
                                   nomatch = length(values) + 1L,
                                   overid = 2L)

      # If .default is scalar, add it to value map

      if (length(.default) <= 1){
        values <- c(values, .default)
        out <- values[key_locs]
      } else {
        out <- values[key_locs]
        no_match_locs <- val_find(key_locs, length(values) + 1L)
        out[no_match_locs] <- .default[no_match_locs]
      }

  } else {
    ## Method with length-n RHS

    # Create a vector filled with `.default` if given, otherwise NA

    if (!is.null(.default)){
      out <- if (length(.default) != 1) .default else rep_len_(.default, N)
    } else {
      if (n_exprs == 0){
        out <- rep_len_(.x[NA_integer_], N)
      } else {
        out <- rep_len_(rhs_list[[1L]][NA_integer_], N)
      }
    }

    for (i in seq_along(exprs)){
      lhs <- lhs_list[[i]]
      rhs <- rhs_list[[i]]
      val_locs <- val_find(.x, lhs)
      out[val_locs] <- if (length(rhs) == 1) rhs else rhs[val_locs]
    }
  }
  out
}
