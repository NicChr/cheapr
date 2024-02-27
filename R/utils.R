is_integerable <- function(x){
  abs(x) <= .Machine$integer.max
}
all_integerable <- function(x, shift = 0){
  all(
    (abs(collapse::frange(x, na.rm = TRUE)) + shift ) <= .Machine$integer.max,
    na.rm = TRUE
  )
}
fill_with_na <- function(x, n = NULL, prop = NULL){
  if (!is.null(n) && !is.null(prop)) {
    stop("either n or prop must be supplied")
  }
  if (!is.null(n)) {
    x[sample.int(length(x), size = n, replace = FALSE)] <- NA
  }
  if (!is.null(prop)) {
    x[sample.int(length(x), size = floor(prop * length(x)),
                 replace = FALSE)] <- NA
  }
  x
}
allv2 <- function(x, value){
  if (!length(x)) {
    return(FALSE)
  }
  collapse::allv(x, value)
}
any_na <- function(x){
  num_na(x) > 0
}
all_na <- function(x){
  num_na(x) == cpp_vector_size(x)
}
