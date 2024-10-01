#' @exportS3Method base::as.double
as.double.integer64 <- function(x, ...){
  cpp_int64_to_double(x)
}
#' @exportS3Method base::as.integer
as.integer.integer64 <- function(x, ...){
  cpp_int64_to_int(x)
}

# collapse methods for integer64
# They are obviously slow but at least result is correct

#' @exportS3Method collapse::fsum
fsum.integer64 <- function(x, g = NULL, ...){

  ## The statement is here because
  ## when there are groups, collapse::fsum
  ## errors when there is a 32-bit int overflow

  if (is.null(g)){
    collapse::fsum(cpp_int64_to_numeric(x), ...)
  } else {
    collapse::fsum(as.double(x), ...)
  }
}
#' @exportS3Method collapse::fmin
fmin.integer64 <- function(x, ...){
  cpp_numeric_to_int64(collapse::fmin(cpp_int64_to_numeric(x), ...))
}
#' @exportS3Method collapse::fmax
fmax.integer64 <- function(x, ...){
  cpp_numeric_to_int64(collapse::fmax(cpp_int64_to_numeric(x), ...))
}
#' @exportS3Method collapse::ffirst
ffirst.integer64 <- function(x, ...){
  cpp_numeric_to_int64(collapse::ffirst(cpp_int64_to_numeric(x), ...))
}
#' @exportS3Method collapse::flast
flast.integer64 <- function(x, ...){
  cpp_numeric_to_int64(collapse::flast(cpp_int64_to_numeric(x), ...))
}
#' @exportS3Method collapse::fmean
fmean.integer64 <- function(x, ...){
  collapse::fmean(cpp_int64_to_numeric(x), ...)
}
#' @exportS3Method collapse::fmedian
fmedian.integer64 <- function(x, ...){
  collapse::fmedian(cpp_int64_to_numeric(x), ...)
}
#' @exportS3Method collapse::fvar
fvar.integer64 <- function(x, ...){
  collapse::fvar(as.double(x), ...)
}
#' @exportS3Method collapse::fsd
fsd.integer64 <- function(x, ...){
  collapse::fsd(as.double(x), ...)
}
#' @exportS3Method collapse::fnth
fnth.integer64 <- function(x, ...){
  collapse::fnth(cpp_int64_to_numeric(x), ...)
}
#' @exportS3Method collapse::fnobs
fnobs.integer64 <- function(x, ...){
  collapse::fnobs(cpp_int64_to_numeric(x), ...)
}
