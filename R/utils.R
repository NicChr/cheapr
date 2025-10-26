#' @noRd

# Like deparse1 but has a cutoff in case of massive strings

deparse2 <- function(expr, collapse = " ", width.cutoff = 500L, nlines = 10L, ...){
  paste(deparse(expr, width.cutoff, nlines = nlines, ...), collapse = collapse)
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

allv2 <- function(x, value){
  if (!length(x)) {
    return(FALSE)
  }
  collapse::allv(x, value)
}

check_length <- function(x, n){
  if (length(x) != n){
    stop(paste(deparse2(substitute(x)), "must have length", n))
  }
}
check_is_df <- function(x){
  if (!inherits(x, "data.frame")){
    stop(paste(deparse2(substitute(x)), "must be a data frame."))
  }
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

list_is_df_like <- function(x){
  collapse::fnunique(list_lengths(x)) <= 1
}
posixlt_is_balanced <- function(x){
  isTRUE(attr(x, "balanced")) && list_is_df_like(x)
}
fill_posixlt <- function(x, classed = TRUE){
  if (!inherits(x, "POSIXlt")){
    stop("x must be a POSIXlt")
  }
  out <- recycle(.args = unclass(x))
  attributes(out) <- attributes(x)
  if (!posixlt_is_balanced(x)){
    attr(out, "balanced") <- NA
  }
  if (!classed){
    class(out) <- NULL
  }
  out
}

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

# The breaks the `cut(x, n)` produces
r_cut_breaks <- function(x, n){
  check_length(n, 1)
  if (is.na(n) || n < 2){
    stop(paste("number of breaks must be >= 2, not", n))
  }
  breaks <- get_breaks(x, n, pretty = FALSE, expand_min = FALSE, expand_max = FALSE)
  nb <- length(breaks)
  adj <- (breaks[nb] - breaks[1]) * 0.001
  breaks[1] <- breaks[1] - adj
  breaks[nb] <- breaks[nb] + adj
  breaks
}

# If args is a plain list of items then extract the first element of
# the top list

as_list_of <- function(...){
  dots <- list(...)
  if (length(dots) == 1 && !is.object(dots[[1L]]) && is.list(dots[[1L]])){
    dots[[1L]]
  } else {
    dots
  }
}

# unevaluated expression list
exprs <- function(...){
  as.list(substitute(alist(...)))[-1L]
}

`%||%` <- function(x, y) if (is.null(x)) y else x
# Sort of the inverse of %||%
`%!||%` <- function(x, y) if (x) NULL else y

# Both below to be safely used in C++ code
fast_match <- function(x, table, nomatch = NA_integer_){
  collapse::fmatch(x, table, overid = 2L, nomatch = nomatch)
}
fast_unique <- function(x){
  unique_(x)
}

vec_setdiff <- function(x, y, unique = FALSE){
  .Call(`_cheapr_cpp_setdiff`, x, y, unique)
}
vec_intersect <- function(x, y, unique = FALSE){
  .Call(`_cheapr_cpp_intersect`, x, y, unique)
}

# `as.numeric` but keep integers as integers
as_numeric <- function(x){
  switch(typeof(x),
         integer = as.integer(x),
         as.double(x))
}

numeric_subtraction <- function(x, y){
  as_numeric(x) - as_numeric(y)
}

numeric_addition <- function(x, y){
  as_numeric(x) + as_numeric(y)
}
