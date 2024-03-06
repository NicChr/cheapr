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
list_as_df <- function(x){
  out <- unclass(x)
  if (length(out) == 0) {
    N <- 0L
  } else {
    N <- length(out[[1L]])
  }
  attr(out, "row.names") <- .set_row_names(N)
  class(out) <- "data.frame"
  out
}
as.character.vctrs_rcrd <- function(x, ...){
  format(x, ...)
}
# Does rcrd fields have rcrds in them?
is_nested_rcrd <- function(x){
  out <- FALSE
  for (i in seq_along(unclass(x))){
    if (inherits(.subset2(x, i), "vctrs_rcrd")){
      out <- TRUE
      break
    }
  }
  out
}
funique.vctrs_rcrd <- function(x, sort = FALSE, ...){
  out <- unique(x, ...)
  if (sort){
    out <- sort(out)
  }
  out
}
funique.POSIXlt <- function(x, ...){
  out <- collapse::funique(list_as_df(x), ...)
  out <- unclass(out)
  attributes(out) <- attributes(x)
  out
}
# funique.vctrs_rcrd <- function(x, sort = FALSE, method = "auto",
#                                decreasing = FALSE, na.last = TRUE, ...){
#   # if (is_nested_rcrd(x)){
#   #  return(unique(x, ...))
#   # }
#   cl <- class(x)
#   out <- collapse::funique(unclass(x), sort = sort, method = method,
#                            decreasing = decreasing, na.last = na.last)
#   class(out) <- cl
#   out
# }
# safe_unique <- function(x, ...){
#   out <- tryCatch(collapse::funique(x, ...), error = function(e) return(".r.error"))
#   if (length(out) == 1 && out == ".r.error"){
#     out <- unique(x, ...)
#   }
#   out
# }
# any_na <- function(x){
#   num_na(x) > 0
# }
# all_na <- function(x){
#   num_na(x) == unlisted_length(x)
# }
