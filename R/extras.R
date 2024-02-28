#' Extra utilities
#'
#' @param x A vector or data frame.
#' @param y A vector or data frame.
#' @param dups Should duplicates be kept? Default is `TRUE`.
#' @param name The column name to assign the names of a vector.
#' @param value The column name to assign the values of a vector.
#' @param breaks See `?cut`.
#' @param labels See `?cut`.
#' @param include.lowest See `?cut`.
#' @param right See `?cut`.
#' @param dig.lab See `?cut`.
#' @param ordered_result See `?cut`.
#' @param table See `?collapse::fmatch`
#' @param ... See `?cut`.
#'
#' @details
#' `intersect_()` and `setdiff_()` are faster and more efficient
#' alternatives to `intersect()` and `setdiff()` respectively. \cr
#' `enframe_()` and `deframe_()` are faster alternatives to
#' `tibble::enframe()` and `tibble::deframe()` respectively. \cr
#' `cut_numeric()` is a faster and more efficient alternative to
#' `cut.default()`.
#'
#' @export
#' @rdname extras
enframe_ <- function(x, name = "name", value = "value"){
  if (inherits(x, "data.frame")) {
    x <- unclass(x)
    attr(x, "row.names") <- NULL
  }
  x_nms <- names(x)
  x <- unname(x)
  if (is.null(x_nms)) {
    out <- list(x)
    names(out) <- value
  }
  else {
    out <- list(x_nms, x)
    names(out) <- c(name, value)
  }
  attr(out, "class") <- c("tbl_df", "tbl", "data.frame")
  attr(out, "row.names") <- .set_row_names(length(x))
  out
}
#' @export
#' @rdname extras
deframe_ <- function(x){
  ncol <- length(names(x))
  if (!(inherits(x, "data.frame") && ncol %in% (1:2))) {
    stop("`x` must be a 1 or 2 col data frame")
  }
  out <- .subset2(x, ncol)
  if (ncol == 2) {
    names(out) <- as.character(.subset2(x, 1L))
  }
  out
}
#' @export
#' @rdname extras
setdiff_ <- function(x, y, dups = TRUE){
  if (!dups){
    x <- collapse::funique(x)
  }
  i <- which_na(collapse::fmatch(x, y, overid = 2L))
  if (inherits(x, "data.frame")){
    x[i, ]
  } else {
    x[i]
  }
}
#' @export
#' @rdname extras
intersect_ <- function(x, y, dups = TRUE){
  if (!dups){
    x <- collapse::funique(x)
  }
  i <- which_not_na(collapse::fmatch(x, y, overid = 2L))
  if (inherits(x, "data.frame")){
    x[i, ]
  } else {
    x[i]
  }
}
#' @export
#' @rdname extras
cut_numeric <- function(x, breaks, labels = NULL, include.lowest = FALSE,
                        right = TRUE, dig.lab = 3L, ordered_result = FALSE, ...){
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (length(breaks) == 1L) {
    if (is.na(breaks) || breaks < 2L)
      stop("invalid number of intervals")
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- collapse::frange(x, na.rm = TRUE))
    if (isTRUE(dx == 0)) {
      dx <- if (rx[1L] != 0)
        abs(rx[1L])
      else 1
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000,
                        length.out = nb)
    }
    else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] +
                               dx/1000)
    }
  }
  else nb <- length(breaks <- sort.int(as.double(breaks)))
  if (anyDuplicated(breaks))
    stop("'breaks' are not unique")
  codes.only <- FALSE
  if (is.null(labels)) {
    for (dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(0 + breaks, digits = dig, width = 1L)
      if (ok <- all(ch.br[-1L] != ch.br[-nb]))
        break
    }
    labels <- if (ok)
      paste0(if (right)
        "("
        else "[", ch.br[-nb], ",", ch.br[-1L], if (right)
          "]"
        else ")")
    else paste0("Range_", seq_len(nb - 1L))
    if (ok && include.lowest) {
      if (right)
        substr(labels[1L], 1L, 1L) <- "["
      else substring(labels[nb - 1L], nchar(labels[nb -
                                                     1L], "c")) <- "]"
    }
  }
  else if (is.logical(labels) && !labels)
    codes.only <- TRUE
  else if (length(labels) != nb - 1L)
    stop("number of intervals and length of 'labels' differ")
  code <- .bincode(x, breaks, right, include.lowest)
  if (!codes.only) {
    levels(code) <- as.character(labels)
    class(code) <- c(if (ordered_result) "ordered" else character(0), "factor")
  }
  code
}
#' @export
#' @rdname extras
`%in_%` <- function(x, table){
  collapse::fmatch(x, table, overid = 2L, nomatch = 0L) > 0L
  # out <- logical(length(x))
  # out[which_not_na(collapse::fmatch(x, table, overid = 2L))] <- TRUE
  # out
}
#' @export
#' @rdname extras
`%!in_%` <- function(x, table){
  is.na(collapse::fmatch(x, table, overid = 2L, nomatch = NA_integer_))
  # out <- logical(length(x))
  # out[which_na(collapse::fmatch(x, table, overid = 2L))] <- TRUE
  # out
}
