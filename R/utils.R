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

# Keep this in-case anyone was using it
fill_with_na <- na_insert

r_cut_breaks <- function(x, n){
  check_length(n, 1)
  stopifnot(n >= 2)
  breaks <- get_breaks(x, n, pretty = FALSE)
  adj <- diff(range(breaks)) * 0.001
  breaks[1] <- breaks[1] - adj
  breaks[length(breaks)] <- breaks[length(breaks)] + adj
  breaks
}

# Is x an atomic type?
# logical, integer, double, character, raw, complex
# including Dates, factors and POSIXcts

is_base_atomic <- function(x){
  (
    is.atomic(x) && (
      !is.object(x) || inherits(x, c("Date", "POSIXct", "factor"))
    )
  ) ||
    is.null(x)
}

combine_factors <- function(...){
  if (nargs() == 0){
    return(NULL)
  }
  if (nargs() == 1){
    dots <- list(...)
    if (!is.object(dots[[1L]]) && is.list(dots[[1L]])){
      return(do.call(combine_factors, dots[[1L]]))
    } else {
      return(dots[[1L]])
    }
  }
  get_levels <- function(x){
    if (is.factor(x)) levels(x) else as.character(x)
  }
  to_char <- function(x){
    if (is.factor(x)) factor_as_character(x) else as.character(x)
  }
  factors <- list(...)
  levels <- lapply(factors, get_levels)

  # Unique combined levels
  new_levels <- collapse::funique(collapse::vec(levels))

  characters <- lapply(factors, to_char)

  # Combine all factor elements (as character vectors)
  factor_(unlist(characters, recursive = FALSE), levels = new_levels,
          na_exclude = !any_na(new_levels))
}

# A very fast 1-D array frequency table
cheapr_table <- function(x, names = TRUE, order = FALSE, na_exclude = FALSE){
  if (is.factor(x)){
    if (na_exclude){
      f <- levels_drop_na(x)
    } else if (any_na(x)){
      f <- levels_add_na(x)
    } else {
      f <- x
    }
  } else {
    f <- factor_(x, order = order, na_exclude = na_exclude)
  }
  lvls <- attr(f, "levels")
  out <- tabulate(f, nbins = length(lvls))
  if (names){
    names(out) <- lvls
  }
  out
}

# Just a wrapper with a cheaper alternative to `c.factor()`
# cheapr_c <- function(..., .check = TRUE){
#   dots <- list(...)
#   if (.check){
#     for (vec in dots){
#       if (is.factor(vec)){
#         return(combine_factors(dots))
#       }
#       if (is.object(vec)){
#         return(do.call(c, dots))
#       }
#     }
#   }
#   `attributes<-`(collapse::vec(dots), NULL)
# }
