#include "cheapr.h"

[[cpp11::register]]
SEXP cpp_replace(SEXP x, SEXP where, SEXP with, bool in_place, bool quiet){

  // Modify `x` in-place?
  bool maybe_shared = MAYBE_SHARED(x);
  if (in_place && maybe_shared && !quiet){
    Rf_warning("`in_place` is `TRUE` but `x` may be shared by multiple objects");
  }

  bool internal_in_place = in_place || !maybe_shared;

  int32_t NP = 0;

  // Clean where vector
  SHIELD(where = clean_locs(where, x)); ++NP;
  const int* RESTRICT p_where = integer_ptr_ro(where);

  // Cast replacement to type of x
  SHIELD(with = cast_(get_r_type(x), with, x)); ++NP;

  R_xlen_t where_size = vec::length(where);
  R_xlen_t with_size = vec::length(with);

  R_xlen_t xi;
  R_xlen_t withi = 0;

  // Shallow copy lists and deep copy data of vectors
  if (!internal_in_place){
    if (Rf_isVectorList(x)){
      SHIELD(x = vec::shallow_copy(x)); ++NP;
    } else {
      SHIELD(x = cpp_semi_copy(x)); ++NP;
    }
  }

  switch (get_r_type(x)){

  case R_null: {
    break;
  }
  case R_lgl:
  case R_int:
  case R_fct: {

    int* RESTRICT p_x = integer_ptr(x);
    const int* RESTRICT p_with = integer_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      p_x[p_where[i] - 1] = p_with[withi];
    }
    break;
  }
  case R_dbl: {

    double* RESTRICT p_x = real_ptr(x);
    const double* RESTRICT p_with = real_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      p_x[p_where[i] - 1] = p_with[withi];
    }
    break;
  }

  case R_int64: {

    int64_t* RESTRICT p_x = integer64_ptr(x);
    const int64_t* RESTRICT p_with = integer64_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      p_x[p_where[i] - 1] = p_with[withi];
    }
    break;
  }

  case R_chr: {

    const r_string_t *p_with = string_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      set_value<r_string_t>(x, p_where[i] - 1, p_with[withi]);
    }
    break;
  }

  case R_cplx: {

    r_complex_t* p_x = complex_ptr(x);
    const r_complex_t *p_with = complex_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      xi = p_where[i] - 1;
      set_value(p_x, xi, p_with[withi]);
    }
    break;
  }

  case R_raw: {

    const r_byte_t *p_with = raw_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      set_value<r_byte_t>(x, p_where[i] - 1, p_with[withi]);
    }
    break;
  }

  case R_date:
  case R_pxt: {
    SEXP x_cls = SHIELD(get_old_class(x)); ++NP;
    attr::set_old_class(x, r_null);
    replace_in_place(x, where, with, false);
    attr::set_old_class(x, x_cls);
    break;
  }

  case R_list: {
    const SEXP *p_with = list_ptr_ro(with);

    for (R_xlen_t i = 0; i < where_size; recycle_index(withi, with_size), ++i){
      SET_VECTOR_ELT(x, p_where[i] - 1, p_with[withi]);
    }
    break;
  }

  case R_df: {

    const SEXP *p_x = list_ptr_ro(x);
    const SEXP *p_with = list_ptr_ro(with);

    int ncol = Rf_length(x);

    if (ncol != Rf_length(with)){
      YIELD(NP);
      Rf_error("`ncol(x)` must equal `ncol(with)`");
    }

    SEXP x_names = SHIELD(get_old_names(x)); ++NP;
    SEXP with_names = SHIELD(get_old_names(with)); ++NP;

    if (!R_compute_identical(x_names, with_names, 0)){
      YIELD(NP);
      Rf_error("Column names must be identical between `x` and `with`");
    }

    for (int j = 0; j < ncol; ++j){
      SET_VECTOR_ELT(x, j, cpp_replace(p_x[j], where, p_with[j], in_place, quiet));
    }

    break;
  }

  case R_unk: {
    SHIELD(x = eval_pkg_fun("base_assign_at", "cheapr", env::base_env, x, where, with)); ++NP;
    break;
  }

  default: {
    YIELD(NP);
    Rf_error("%s cannot handle an object of R type %s", __func__, r_type_char(x));
  }
  }
  YIELD(NP);
  return x;
}

void replace_in_place(SEXP x, SEXP where, SEXP with, bool quiet){
  cpp_replace(x, where, with, true, quiet);
}
