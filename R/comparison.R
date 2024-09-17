# cheapr_equals <- function(x, y){
#   UseMethod("cheapr_equals")
# }
# cheapr_equals.default <- function(x, y){
#   x == y
# }
# cheapr_equals.data.frame <- function(x, y){
#   cpp_check_nested_lengths(x, y)
#   cpp_data_frames_equal(x, y)
# }

# cheapr_equals <- function(x, y){
#   if (is.list(x) && is.list(y)){
#     recycled <- recycle(x, y)
#     x <- recycled[[1]]
#     y <- recycled[[2]]
#     cpp_check_nested_lengths(x, y)
#     common <- cpp_cast_common(x, y)
#     cpp_data_frames_equal2(common[[1]], common[[2]])
#   } else {
#     x == y
#   }
# }


# cheapr_cast <- function(x, y, .keep_attrs = FALSE){
#   if (identical(class(x), class(y))){
#     x
#   } else {
#     x[0] <- y[0]
#     if (.keep_attrs){
#       `mostattributes<-`(x, attributes(y))
#     }
#     x
#   }
# }

character_compare <- function(x, y, .op){
  if (!( is.character(x) || is.character(y) )){
    stop("Either x or y  must be a character vector")
  }
  funs <- c("==", ">", "<", ">=", "<=", "!=", "=")
  fun <- if (is.character(.op)) .op else as.character(substitute(.op))
  type <- match(fun, funs)
  op <- switch(funs[type],
               `==` = `==`,
               `>` = `>`,
               `<` = `<`,
               `>=` = `>=`,
               `<=` = `<=`,
               `!=` = `!=`,
               `=` = stop("`=` cannot be used for comparisons, perhaps you meant `==`?"),
               stop("Unsupported comparison"))
  cpp_character_compare(x, y, type)
}

# Comparison operators for character vectors
#
# Currently (and this is unlikely to change), only `==` and `!=` is supported
# as cheapr is unicode un-aware and therefore can't make comparisons in
# different locales. For this type of comparison, the 'stringi' package
# can be used.
# These operators are somewhat cheaper when `x` or `y` is numeric and the
# other is a character vector because the entire vector isn't coerced
# to character. For floating-point vectors it's likely the outputs
# will not match base R operator outputs.
#
#
# A logical vector reflecting the requested comparison.
#
`%ch==%` <- function(x, y){
  cpp_character_compare(x, y, 1L)
}
`%ch!=%` <- function(x, y){
  cpp_character_compare(x, y, 6L)
}
`%ch>%` <- function(x, y){
  character_compare(x, y, ">")
}
`%ch<%` <- function(x, y){
  character_compare(x, y, "<")
}
`%ch>=%` <- function(x, y){
  character_compare(x, y, ">=")
}
`%ch<=%` <- function(x, y){
  character_compare(x, y, "<=")
}
