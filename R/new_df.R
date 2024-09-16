#' Fast data frame constructor
#'
#' @param ... Key-value pairs.
#' @param .nrows `integer(1)` (Optional) number of rows. \cr
#' Commonly used to initialise a 0-column data frame with rows.
#' @param .recycle `logical(1)` Should arguments be recycled?
#' Default is `FALSE`.
#' @param .name_repair `logical(1)` Should duplicate names be made unique?
#' Default is `FALSE`.
#'
#' @returns A `data.frame`
#'
#' @export
new_df <- function(..., .nrows = NULL, .recycle = FALSE, .name_repair = FALSE){

  out <- list_named(...)

  # Recycle
  if (.recycle){
    out <- do.call(function(...) recycle(..., length = .nrows), out)
  }

  if (is.null(.nrows)){
    if (length(out) == 0L){
      row_names <- integer()
    } else {
      N <- NROW(.subset2(out, 1L))
      row_names <- c(NA_integer_, -N)
    }
  } else {
    row_names <- .set_row_names(.nrows)
  }

  out_names <- as.character(attr(out, "names", TRUE))

  if (.name_repair){
    out_names <- unique_name_repair(out_names)
  }

  attr(out, "names") <- out_names
  attr(out, "row.names") <- row_names
  class(out) <- "data.frame"
  out
}
unique_name_repair <- function(x, .sep = "..."){
  if (is.null(x)) {
    return(x)
  }
  x <- as.character(x)
  col_seq <- seq_along(x)
  which_dup <- which_(collapse::fduplicated(x, all = TRUE))
  x[which_dup] <- paste0(x[which_dup], .sep, col_seq[which_dup])
  x
}
