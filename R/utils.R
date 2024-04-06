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

list_as_df <- cpp_list_as_df

df_as_tbl <- function(x){
  out <- list_as_df(x)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
as.character.vctrs_rcrd <- function(x, ...){
  format(x, ...)
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
check_length <- function(x, n){
  if (length(x) != n){
    stop(paste(deparse1(substitute(x)), "must have length", n))
  }
}
check_is_df <- function(x){
  if (!inherits(x, "data.frame")){
    stop(paste(deparse1(substitute(x)), "must be a data frame."))
  }
}
df_add_cols <- function(data, cols){
  nms <- names(cols)
  if (is.null(nms)){
    stop("cols must be a named list")
  }
  for (i in seq_along(cols)){
    data[[nms[i]]] <- cols[[i]]
  }
  data
}
which_in <- function(x, table){
  which_not_na(collapse::fmatch(x, table, overid = 2L, nomatch = NA_integer_))
}
which_not_in <- function(x, table){
  which_na(collapse::fmatch(x, table, overid = 2L, nomatch = NA_integer_))
}
tzone <- function(x){
  out <- attr(x, "tzone")
  if (is.null(out)) {
    ""
  }
  else {
    out[[1]]
  }
}

# Recycle arguments
recycle <- function (..., length = NULL){
  out <- cpp_list_rm_null(list(...))
  lens <- lengths_(out)
  uniq_lens <- collapse::fnunique(lens)
  if (is.null(length)) {
    if (length(lens)) {
      N <- max(lens)
    }
    else {
      N <- 0L
    }
  }
  else {
    N <- length
  }
  N <- N * (!collapse::anyv(lens, 0L))
  recycle <- which_(lens != N)
  out[recycle] <- lapply(out[recycle], rep_len, N)
  out
}
n_dots <- function(...){
  nargs()
}
set_attr <- cpp_set_add_attr
set_attrs <- cpp_set_add_attributes
set_rm_attr <- cpp_set_rm_attr
set_rm_attrs <- cpp_set_rm_attributes

balance_posixlt <- function(x){
  balance_pos <- tryCatch(get("balancePOSIXlt",
                              asNamespace("base"),
                              inherits = FALSE),
                          error = function(e) return(".r.error"))
  if (is.character(balance_pos) && length(balance_pos) == 1 && balance_pos == ".r.error"){
    unclass(x)
  } else {
    balance_pos(x, fill.only = FALSE, classed = FALSE)
  }
}
