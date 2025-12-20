threads_are_valid <- function(n){
  if (!is.numeric(n)){
    return(FALSE)
  }
  if (length(n) != 1){
    return(FALSE)
  }

  if (!isTRUE(n >= 1)){
    return(FALSE)
  }
  TRUE
}

#' Get and set the number of OpenMP threads to be used in cheapr
#'
#' @description
#' The default number of threads in
#' cheapr is 2 (if your system has at least 2 available threads).
#' You can change the number of threads via `set_threads()`.
#' To see the number of threads currently being used, run `get_threads()`.
#' To see the maximum number of threads your machine can use,
#' run `get_max_threads()`
#'
#' @param n `[integer(1)]` - Number of threads to use.
#'
#'
#' @returns
#' `get_max_threads()` returns the max number of threads you can use. \cr
#' `get_threads()` returns the number of threads that are being used. \cr
#' `set_threads()` invisibly sets the number of threads to be used.
#'
#' @rdname threads
#' @export
get_max_threads <- cpp_max_threads
#' @rdname threads
#' @export
set_threads <- function(n){
  if (!threads_are_valid(n)){
    stop(paste("Please supply a valid number of threads between 1 and", get_max_threads()))
  }
  threads <- min(as.integer(n), get_max_threads())
  threads <- max(1L, threads)
  options(cheapr.cores = threads)
}
#' @rdname threads
#' @export
get_threads <- function(){
  threads <- getOption("cheapr.cores", default = 2L)
  if (!threads_are_valid(threads)){
    warning("Invalid threads detected, setting threads to default 2")
    threads <- 2L
    set_threads(threads)
  }
  threads
}

