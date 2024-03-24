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
df_select <- function(x, i){
  attrs <- attributes(x)
  out <- cpp_list_rm_null(unclass(x)[i])
  attrs[["names"]] <- attr(out, "names")
  attrs[["row.names"]] <- .row_names_info(x, type = 0L)
  attributes(out) <- attrs
  out
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

# Efficient data frame subset
# With the exception of which_() this is surprisingly all base R...
# It relies on sset which falls back on `[` when no method is found.
df_subset <- function(x, i = 0, j = seq_along(x)){
  nrows <- length(attr(x, "row.names"))
  if (is.logical(i)){
    check_length(i, nrows)
    i <- which_(i)
  }
  l <- as.list(
    df_select(x, j)
  )
  if (length(l) == 0){
    out <- list_as_df(l)
    attr(out, "row.names") <- .set_row_names(
      length(
        seq_len(nrows)[i]
      )
    )
  } else {
    out <- list_as_df(
      lapply(l, sset, i)
    )
  }
  out
}

# safe_unique <- function(x, ...){
#   out <- tryCatch(collapse::funique(x, ...), error = function(e) return(".r.error"))
#   if (length(out) == 1 && out == ".r.error"){
#     out <- unique(x, ...)
#   }
#   out
# }
