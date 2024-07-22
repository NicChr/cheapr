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
