#' Utilities for creating many sequences
#'
#' @description
#' `sequence_` is an extension to [sequence] which
#' accepts decimal number increments. \cr
#' `seq_id` can be paired with `sequence_` to group individual sequences. \cr
#' `seq_` is a vectorised version of [seq]. \cr
#' `window_sequence` creates a vector of window sizes for rolling calculations. \cr
#' `lag_sequence` creates a vector of lags for rolling calculations. \cr
#' `lead_sequence` creates a vector of leads for rolling calculations. \cr
#'
#' @param size Vector of sequence lengths.
#' @param from Start of sequence(s).
#' @param to End of sequence(s).
#' @param add_id Should the ID numbers of the sequences be added as names?
#' Default is `FALSE`.
#' @param as_list Should a list of sequences be returned?
#' Setting to `TRUE` would
#' place each distinct sequence vector into a distinct list element.
#' The default is `FALSE`.
#'
#' @param by Unit increment of sequence(s).
#' @param k Window/lag size.
#' @param partial Should partial windows/lags be returned? Default is `TRUE`.
#' @param ascending Should window sequence be ascending? Default is `TRUE`.
#'
#' @returns
#' A vector of length `sum(size)` except for `seq_` which
#' returns a vector of size `sum((to - from) / (by + 1))`
#'
#' @details
#' `sequence_()` works in the same way as `sequence()` but can accept
#' non-integer `by` values.
#' It also recycles `from` and `to`, in the same way as `sequence()`. \cr
#' If any of the sequences contain values > `.Machine$integer.max`,
#' then the result will always be a double vector.
#'
#' `from` can be also be a date, date-time, or any object that supports
#' addition and multiplication.
#'
#' `seq_()` is a vectorised version of `seq()` that strictly accepts
#' only the arguments `from`, `to` and `by`. \cr
#'
#' @examples
#' library(cheapr)
#'
#' # These two functions are similar
#' sequence(1:3);sequence_(1:3)
#'
#' # sequence_() can handle any numeric vector sequence
#' sequence(1:3, by = 0.1);sequence_(1:3, by = 0.1)
#'
#' # Alternatively return as a list of sequences
#'
#' sequence_(1:3, by = 0.1, as_list = TRUE)
#'
#' # Add IDs to the sequences
#' sequence_(1:3, by = 0.1, add_id = TRUE)
#' # Turn this quickly into a data frame
#' seqs <- sequence_(1:3, by = 0.1, add_id = TRUE)
#' new_df(name = names(seqs), seq = seqs)
#'
#' sequence(c(3, 2), by = c(-0.1, 0.1));sequence_(c(3, 2), by = c(-0.1, 0.1))
#'
#' # Vectorised version of seq()
#' seq_(1, 10, by = c(1, 0.5))
#' # Same as above
#' c(seq(1, 10, 1), seq(1, 10, 0.5))
#'
#' # Again, as a list of sequences
#' # 2 different start points and 2 different increments
#' seq_(from = c(-1, 1), 3, by = c(1, 0.5), as_list = TRUE)
#'
#' # Programmers may use seq_size() to determine final sequence lengths
#'
#' sizes <- seq_size(1, 10, by = c(1, 0.5))
#' print(paste(c("sequence sizes: (", sizes, ") total size:", sum(sizes)),
#'             collapse = " "))
#'
#' # Or return as a list of sequences
#' # Note that these lengths will match the above line of code
#' seq_(1, 10, by = c(1, 0.5), as_list = TRUE) |>
#'   list_lengths()
#'
#' # Sequences of dates with different increments
#'
#' from <- Sys.Date()
#' to <- from + 10
#' by <- c(1, 2, 3)
#' date_seqs <- seq_(from, to, by, as_list = TRUE)
#' lapply(date_seqs, \(x) `class<-`(x, "Date"))
#'
#' # Utilities for rolling calculations
#'
#' # A window sequence of size 3 for a vector of size 10
#' # This tells us how big the window should be when looking backwards
#' window_sequence(10, 3, partial = FALSE)
#' window_sequence(10, 3, partial = TRUE)
#'
#' window_sequence(c(3, 5), 3)
#' window_sequence(c(3, 5), 3, partial = FALSE)
#' window_sequence(c(3, 5), 3, partial = TRUE, ascending = FALSE)
#'
#' # Lag sequence of size 3 for a vector of size 10
#' # This tells us how for we should look backwards at any given point
#' lag_sequence(10, 3, partial = FALSE)
#' # How far to look forwards
#' lead_sequence(10, 3, partial = FALSE)
#'
#' lag_sequence(10, 3, partial = TRUE)
#' lead_sequence(10, 3, partial = TRUE)

#' # One can for example use these in data.table::frollsum
#'
#' @rdname sequences
#' @export
sequence_ <- function(size, from = 1L, by = 1L, add_id = FALSE, as_list = FALSE){
  out <- cpp_sequence(as.integer(size), from, by, as.logical(as_list))
  if (add_id){
    if (as_list){
      names(out) <- seq_along(out)
    } else {
      names(out) <- seq_id(size)
    }
  }
  out
}
#' @rdname sequences
#' @export
seq_id <- cpp_sequence_id
# Alternative base R way of defining seq_id
# rep.int(seq_along(size), times = size)

#' @rdname sequences
#' @export
seq_ <- function(from = 1L, to = 1L, by = 1L, add_id = FALSE, as_list = FALSE){
  out_size <- seq_size(from = from, to = to, by = by)
  sequence_(out_size, from = from, by = by, add_id = add_id, as_list = as_list)
}
#' @rdname sequences
#' @export
seq_size <- function(from, to, by = 1L){
  del <- unclass(to - from)
  if (is.integer(by) && allv2(by, 1L)){
    size <- del
  } else {
    size <- del / by
    size[which_val(del, 0)] <- 0
  }
  size_rng <- collapse::frange(size, na.rm = TRUE)
  if (isTRUE(any(size_rng < 0))){
    stop("At least 1 sequence length is negative, please check the sign of by")
  }
  if (length(size) == 0 || all(is_integerable(abs(size_rng) + 1), na.rm = TRUE)){
    if (is.integer(size)){
      size <- set_add(size, 1L)
    } else {
      # size <- as.integer(size + 1e-10) + 1L
      size <- set_add(size, 1e-10)
      size <- as.integer(size)
      size <- set_add(size, 1L)
    }
  } else {
    size <- trunc(size + 1e-10) + 1
  }
  size
}
seq_from <- function(size, to, by = 1L){
  to - (pmax(size - 1L, 0L) * by)
}
seq_to <- function(size, from, by = 1L){
  from + (pmax(size - 1L, 0L) * by)
}
seq_by <- function(size, from, to){
  del <- as.double((to - from))
  out <- del / pmax.int(size - 1L, 0L)
  out[which_val(del, 0)] <- 0
  out
}
#' @rdname sequences
#' @export
window_sequence <- function(size, k, partial = TRUE, ascending = TRUE, add_id = FALSE) {
  out <- cpp_window_sequence(size, k, partial, ascending)
  if (add_id){
    names(out) <- seq_id(size)
  }
  out
}
#' @rdname sequences
#' @export
lag_sequence <- function(size, k, partial = TRUE, add_id = FALSE) {
  out <- cpp_lag_sequence(size, k, partial)
  if (add_id){
    names(out) <- seq_id(size)
  }
  out
}
#' @rdname sequences
#' @export
lead_sequence <- function(size, k, partial = TRUE, add_id = FALSE) {
  out <- cpp_lead_sequence(size, k, partial)
  if (add_id){
    names(out) <- seq_id(size)
  }
  out
}
