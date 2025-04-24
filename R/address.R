#' Memory address of R object
#'
#' @param x An R object.
#'
#' @returns
#' Memory address of R object.
#'
#' @export
address <- function(x){
  .Call(`_cheapr_cpp_address`, x)
}
