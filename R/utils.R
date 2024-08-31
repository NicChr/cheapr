#' @noRd

# Include this as it's not available in older versions of R
deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...){
  paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

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

cheapr_rep_len <- function(x, length.out){
  if (inherits(x, "data.frame")){
    sset(x, rep_len(attr(x, "row.names"), length.out))
  } else {
    rep(x, length.out = length.out)
  }
}

set_attr <- cpp_set_add_attr
set_attrs <- cpp_set_add_attributes
set_rm_attr <- cpp_set_rm_attr
set_rm_attrs <- cpp_set_rm_attributes

cpp_list_rm_null <- function(x, always_shallow_copy = TRUE){
  cpp_drop_null(x, always_shallow_copy)
}
list_is_df_like <- function(x){
  collapse::fnunique(lengths_(x)) <= 1
}
posixlt_is_balanced <- function(x){
  isTRUE(attr(x, "balanced")) && list_is_df_like(x)
}
fill_posixlt <- function(x, classed = TRUE){
  if (!inherits(x, "POSIXlt")){
    stop("x must be a POSIXlt")
  }
  out <- unclass(x)
  out <- do.call(recycle, out)
  attributes(out) <- attributes(x)
  if (!posixlt_is_balanced(x)){
    attr(out, "balanced") <- NA
  }
  if (!classed){
    class(out) <- NULL
  }
  out
}

# balance_posixlt <- function(x, fill.only = FALSE, classed = TRUE){
#   balance_pos <- tryCatch(get("balancePOSIXlt",
#                               asNamespace("base"),
#                               inherits = FALSE),
#                           error = function(e) return(".r.error"))
#   if (is.character(balance_pos) && length(balance_pos) == 1 && balance_pos == ".r.error"){
#     fill_posixlt(x, classed = classed)
#   } else {
#     balance_pos(x, fill.only = fill.only, classed = classed)
#   }
# }

#' @exportS3Method base::as.character
as.character.vctrs_rcrd <- function(x, ...){
  format(x, ...)
}
#' @exportS3Method collapse::funique
funique.vctrs_rcrd <- function(x, sort = FALSE, ...){
  out <- unique(x, ...)
  if (sort){
    out <- sort(out)
  }
  out
}
#' @exportS3Method collapse::funique
funique.POSIXlt <- function(x, sort = FALSE, ...){
  out <- fill_posixlt(x, classed = FALSE)
  # out <- balance_posixlt(x, fill.only = TRUE, classed = FALSE)
  out_attrs <- attributes(out)
  out <- list_as_df(out)
  out <- collapse::funique(out, sort = FALSE)
  if (sort){
    o <- collapse::radixorderv(
      sset(out, j = c("year", "yday", "hour", "min", "sec")),
      ...
    )
    out <- sset(out, o)
  }
  attributes(out) <- out_attrs
  class(out) <- class(x)
  out
}

n_dots <- function(...){
  nargs()
}

dot_expr_names <- function(...){
  vapply(substitute(alist(...))[-1L], deparse1, "", USE.NAMES = FALSE)
}

named_dots <- function(...){
  dots <- list(...)

  dot_nms <- names(dots)

  if (is.null(dot_nms)){
    names(dots) <- dot_expr_names(...)
  } else if (!all(nzchar(dot_nms))){
    empty <- which_(!nzchar(dot_nms))
    expr_names <- dot_expr_names(...)
    dot_nms[empty] <- expr_names[empty]
    names(dots) <- dot_nms
  }
  dots
}
