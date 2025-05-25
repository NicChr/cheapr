
# rebuild.R has another dependency

set_attr <- cpp_set_add_attr
set_attrs <- cpp_set_add_attributes
set_rm_attr <- cpp_set_rm_attr
set_rm_attrs <- cpp_set_rm_attributes

cpp_list_rm_null <- function(x, always_shallow_copy = TRUE){
  cpp_drop_null(x, always_shallow_copy)
}

# To not break dependency
cpp_reconstruct <- cpp_rebuild

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

# Keep this in-case anyone was using it
fill_with_na <- na_insert

# Keep this for fastplyr otherwise it breaks dependency
df_select <- function(x, j = NULL){
  cpp_rebuild(
    sset_col(x, missing(j) %!||% j), x,
    "names", val_rm(names(attributes(x)), "names"),
    TRUE
  )
}

# Kept for reverse compatibility reasons
cpp_sset_df <- sset_row
