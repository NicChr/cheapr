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

## Like the R equivalents
## except an attr '.n.true' is returned
## showing the number of true values in lgl vec
## This is very useful for further subsetting
## When '.n.true' is 0, subsetting is very fast compare to base


cheapr_compare <- function(x, y, .op){
  if (!is.object(x) && is.atomic(x) && !is.object(y) && is.atomic(y)){
    cpp_compare(
      x, y,
      switch(if (is.character(.op)) .op else as.character(substitute(.op)),
             `==` = 1L,
             `>` = 2L,
             `<` = 3L,
             `>=` = 4L,
             `<=` = 5L,
             `!=` = 6L,
             0L)
    )
  } else {
    match.fun(.op)(x, y)
  }
}

`%v==%` <- function(e1, e2){
  cheapr_compare(e1, e2, "==")
}
`%v!=%` <- function(e1, e2){
  cheapr_compare(e1, e2, "!=")
}
`%v>%` <- function(e1, e2){
  cheapr_compare(e1, e2, ">")
}
`%v<%` <- function(e1, e2){
  cheapr_compare(e1, e2, "<")
}
`%v>=%` <- function(e1, e2){
  cheapr_compare(e1, e2, ">=")
}
`%v<=%` <- function(e1, e2){
  cheapr_compare(e1, e2, "<=")
}
