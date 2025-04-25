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

list_as_df <- cpp_list_as_df

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

# Keep this in-case anyone was using it
fill_with_na <- na_insert

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

# Combine levels of factors
# Converts non factors into character vectors
combine_levels <- function(...){
  cpp_combine_levels(list(...))
}

combine_factors <- function(...){
  cpp_combine_factors(list(...))
}


# Keeping this as other packages may use it
df_add_cols <- function(data, cols){
  if (!(is.list(cols) && !is.null(names(cols)))){
    stop("cols must be a named list")
  }
  N <- length(attr(data, "row.names"))
  out <- unclass(data)
  temp <- unclass(cols)
  for (col in names(temp)){
    out[[col]] <- if (is.null(temp[[col]])) NULL else cheapr_rep_len(temp[[col]], N)
  }
  class(out) <- class(data)
  out
}

# unevaluated expression list
exprs <- function(...){
  as.list(substitute(alist(...)))[-1]
}

# Sort of the inverse of %||%
`%!||%` <- function(x, y) if (x) NULL else y

# Turn negative indices to positives
neg_indices_to_pos <- function(exclude, n){
  if (n == 0){
    integer()
  } else {
    which_not_in(
      seq.int(from = -1L, to = -as.integer(n), by = -1L),
      as.integer(exclude)
    )
  }
}

# Both below to be safely used in C++ code
fast_match <- function(x, table, nomatch = NA_integer_){
  collapse::fmatch(x, table, overid = 2L, nomatch = nomatch)
}
fast_unique <- function(x){
  collapse::funique(x)
}

vec_setdiff <- function(x, y, unique = FALSE){
  .Call(`_cheapr_cpp_setdiff`, x, y, unique)
}
vec_intersect <- function(x, y, unique = FALSE){
  .Call(`_cheapr_cpp_intersect`, x, y, unique)
}
