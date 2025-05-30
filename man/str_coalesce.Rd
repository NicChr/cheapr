% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/str_coalesce.R
\name{str_coalesce}
\alias{str_coalesce}
\title{Coalesce character vectors}
\usage{
str_coalesce(..., .args = NULL)
}
\arguments{
\item{...}{Character vectors to coalesce.}

\item{.args}{An alternative to \code{...} for easier programming with lists.}
}
\value{
A coalesced character vector of length corresponding to the recycled
size of supplied character vectors. See \code{?recycle} for details.
}
\description{
\code{str_coalesce()} find the first non empty string \code{""}.
This is particularly useful for assigning and fixing the
names of R objects.

In this implementation, the empty string \code{""} has priority over
\code{NA} which means \code{NA} is only returned when all
values are \code{NA}, e.g. \code{str_coalesce(NA, NA)}.
}
\details{
\code{str_coalesce(x, y)} is equivalent to
\code{if_else(x != "" & !is.na(x), x, y)}.
}
\examples{
library(cheapr)

# Normal examples
str_coalesce("", "hello")
str_coalesce("", NA, "goodbye")

# '' always preferred
str_coalesce("", NA)
str_coalesce(NA, "")

# Unless there are only NAs
str_coalesce(NA, NA)

# `str_coalesce` is vectorised

x <- val_insert(letters, "", n = 10)
y <- val_insert(LETTERS, "", n = 10)

str_coalesce(x, y)

# Using `.args` instead of `do.call` is much more efficient
library(bench)
x <- cheapr_rep_len(list(letters), 10^3)

mark(do.call(str_coalesce, x),
     str_coalesce(.args = x),
     iterations = 50)

}
