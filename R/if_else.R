#' Cheaper version of `ifelse()`
#'
#' @param condition [logical] A condition which will be used to
#' evaluate the if else operation.
#' @param true Value(s) to replace `TRUE` instances.
#' @param false Value(s) to replace `FALSE` instances.
#' @param na Catch-all value(s) to replace all other instances,
#' where `is.na(condition)`.
#'
#' @seealso [case] [val_match]
#'
#' @returns
#' A vector the same length as condition,
#' using a common type between `true`, `false` and `default`.
#'
#' @export
cheapr_if_else <- function(condition, true, false, na = false[NA_integer_]){

  if (!is.logical(condition)){
    stop("condition must be a logical vector")
  }
  if (length(true) != 1 && length(true) != length(condition)){
    stop("`length(true)` must be 1 or `length(condition)`")
  }
  if (length(false) != 1 && length(false) != length(condition)){
    stop("`length(false)` must be 1 or `length(condition)`")
  }
  if (length(na) != 1 && length(na) != length(condition)){
    stop("`length(na)` must be 1 or `length(condition)`")
  }

  if (is.factor(true) || is.factor(false) || is.factor(na)){
    template <- combine_factors(true[1L], false[1L], na[1L])[0L]
  } else {
    template <- c(true[1L], false[1L], na[1L])[0L]
  }

  true <- cast(true, template)
  false <- cast(false, template)
  na <- cast(na, template)

  if (is_base_atomic(true) && is_base_atomic(false) && is_base_atomic(na)){
    return(`mostattributes<-`(
      cpp_if_else(condition, true, false, na),
      attributes(template)
    ))
  }

  # Catch-all method

  lgl_val_counts <- cpp_lgl_count(condition)
  n_true <- lgl_val_counts["true"]
  n_false <- lgl_val_counts["false"]
  n_na <- lgl_val_counts["na"]

  if (n_true == length(condition)){
    if (length(true) == 1){
      return(rep(true, length(condition)))
    } else {
      return(true)
    }
  }

  if (n_false == length(condition)){
    if (length(false) == 1){
      return(rep(false, length(condition)))
    } else {
      return(false)
    }
  }

  if (n_na == length(condition)){
    if (length(na) == 1){
      return(rep(na, length(condition)))
    } else {
      return(na)
    }
  }

  # The else part is most likely to be most prominent
  if (length(false) == length(condition)){
    out <- false
  } else {
    out <- rep(false, length.out = length(condition))
  }

  lgl_locs <- cpp_lgl_locs(condition, n_true = n_true, n_false = n_false,
                           include_true = TRUE, include_false = FALSE,
                           include_na = TRUE)
  true_locs <- lgl_locs[["true"]]
  na_locs <- lgl_locs[["na"]]

  if (length(true) == 1){
    out[true_locs] <- true
  } else {
    out[true_locs] <- true[true_locs]
  }
  if (length(na) == 1){
    out[na_locs] <- na
  } else {
    out[na_locs] <- na[na_locs]
  }
  out
}

# cheapr_case_when <- function(..., .default = NULL){
#   # key_value_list <- lapply(substitute(alist(...))[-1], as.character)
#   conditions <- as.list(match.call())[-1L]
#   default_val <- eval(conditions[[".default"]], envir = parent.frame())
#   conditions <- conditions[-match(".default", names(conditions),
#                                   nomatch = length(conditions) + 1L)]
#   n_conditions <- length(conditions)
#
#   if (n_conditions == 0){
#     stop("Please supply >= 1 case when statements")
#   }
#
#   rhs_list <- new_list(n_conditions)
#   lgl_list <- new_list(n_conditions)
#
#   # Evaluate first condition
#
#   condition <- conditions[[1L]]
#   if (as.character(condition[[1L]]) != "~" || length(condition) != 3){
#     stop("Please supply formula expressions")
#   }
#   rhs <- eval(condition[[3L]], envir = parent.frame())
#   rhs_list[[1L]] <- rhs
#
#   ## TO-do CHECK length of rhs
#
#   # template <- rhs[1L]
#
#   lgl_subset <- seq_along(condition)
#   lgl <- eval(condition[[2L]], envir = parent.frame())
#   true_locs <- val_find(lgl, TRUE)
#   if (length(rhs) == 1){
#     out <- rep(rhs[NA_integer_], length.out = length(lgl))
#   } else {
#     # out <- recycle(rhs, length = max(length(rhs), length(lgl)))[[1L]]
#     out <- rhs
#     rhs <- rhs[true_locs]
#   }
#   out[true_locs] <- rhs
#
#   if (n_conditions >= 2){
#     # lgl2 <- lgl
#     lgl_or <- r_copy(lgl)
#     # lgl_and <- r_copy(lgl)
#     lgl3 <- logical(length(lgl))
#     for (i in 2:n_conditions){
#
#       condition <- conditions[[i]]
#       if (as.character(condition[[1L]]) != "~" || length(condition) != 3){
#         stop("Please supply formula expressions")
#       }
#       rhs <- eval(condition[[3L]], envir = parent.frame())
#
#       ## TO-DO: Check length of rhs here
#       lgl <- eval(condition[[2L]], envir = parent.frame())
#
#       cpp_set_copy_elements(source = lgl, target = lgl3)
#       # cpp_set_replace(lgl3, val_find(lgl2, TRUE), FALSE)
#       # lgl2 <- lgl2 | lgl
#
#       cpp_set_replace(lgl3, val_find(lgl_or, TRUE), FALSE)
#       cpp_set_or(lgl_or, lgl)
#       # cpp_set_and(lgl_and, lgl)
#
#
#       true_locs <- val_find(lgl3, TRUE)
#       if (length(true_locs) > 0){
#         if (length(rhs) == 1){
#           out[true_locs] <- rhs
#         } else {
#           out[true_locs] <- rhs[true_locs]
#         }
#       }
#     }
#   }
#
#
#   if (!is.null(default_val)){
#     default_locs <- val_find(lgl_or, TRUE, invert = TRUE)
#
#     if (length(default_locs) > 0){
#       if (length(default_val) == 1){
#         out[default_locs] <- default_val
#       } else {
#         out[default_locs] <- default_val[default_locs]
#       }
#     }
#
#   }
#   out
# }

# cheapr_case_when <- function(..., .default = NA){
#   # key_value_list <- lapply(substitute(alist(...))[-1], as.character)
#   conditions <- as.list(match.call())[-1L]
#   n_conditions <- length(conditions)
#
#   if (n_conditions == 0){
#     stop("Please supply >= 1 case when statements")
#   }
#
#   rhs_list <- new_list(n_conditions)
#
#   # Evaluate first condition
#
#   condition <- conditions[[1L]]
#   if (as.character(condition[[1L]]) != "~" || length(condition) != 3){
#     stop("Please supply formula expressions")
#   }
#   rhs <- eval(condition[[3L]], envir = parent.frame())
#   rhs_list[[1L]] <- rhs
#
#   ## TO-do CHECK length of rhs
#
#   # template <- rhs[1L]
#
#   lgl_subset <- seq_along(condition)
#   lgl <- eval(condition[[2L]], envir = parent.frame())
#   true_locs <- val_find(lgl, TRUE)
#   if (length(rhs) == 1){
#     out <- rep(rhs[NA_integer_], length.out = length(lgl))
#   } else {
#     out <- rhs
#     rhs <- rhs[true_locs]
#   }
#   out[true_locs] <- rhs
#
#   if (n_conditions >= 2){
#     # temp <- out
#     for (i in 2:n_conditions){
#
#       condition <- conditions[[i]]
#       if (as.character(condition[[1L]]) != "~" || length(condition) != 3){
#         stop("Please supply formula expressions")
#       }
#       rhs <- eval(condition[[3L]], envir = parent.frame())
#
#       ## TO-DO: Check length of rhs here
#
#       lgl_subset <- val_find(lgl, TRUE, invert = TRUE)
#       lgl <- eval(condition[[2L]], envir = parent.frame())
#       true_locs <- val_find(lgl[lgl_subset], TRUE)
#       if (length(true_locs) > 0){
#         if (length(rhs) != 1){
#           rhs <- rhs[true_locs]
#         }
#         out[lgl_subset][true_locs] <- rhs
#       }
#     }
#   }
#
#   # The above currently doesn't work
#
#   # The below pseudo code will work
#
#   # lgl1 <- x %in% 1:4
#   # lgl2 <- x %in% 2:3
#   # lgl3 <- x %in% 3:5
#   #
#   # lgl <- lgl2
#   # lgl[lgl1] <- FALSE
#   # # lgl[lgl2] <- TRUE
#   #
#   # lgl <- lgl3
#   # lgl[lgl1 | lgl2] <- FALSE
#   # lgl[lgl3] <- TRUE
#
#   # Basically we have a logical vector template
#   # As the current logical condition
#   # Using previous conditions with an OR operator
#   # We set all those to FALSE
#
#   # if (!(length(.default) == 1 && is.na(.default)) && length(lgl_subset) != 0){
#   #   default_locs <- val_find(lgl[lgl_subset], TRUE, invert = TRUE)
#   #
#   #   if (length(default_locs) > 0){
#   #     if (length(.default) != 1){
#   #
#   #     }
#   #   }
#   #
#   # }
#   out
# }
