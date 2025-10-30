#include "cheapr.h"

[[cpp11::register]]
SEXP cpp_assign(SEXP x, SEXP where, SEXP with, bool in_place){

  // Modify `x` in-place?
  bool maybe_shared = MAYBE_SHARED(x);
  if (in_place && maybe_shared){
    Rf_warning("`in_place` is `TRUE` but `x` may be shared by multiple objects");
  }

  bool internal_in_place = in_place || !maybe_shared;

  int32_t NP = 0;

  // Convere where to integer
  SHIELD(where = cast<r_integer_t>(where, R_NilValue)); ++NP;
  const int* RESTRICT p_where = INTEGER_RO(where);

  // Cast replacement to type of x
  SHIELD(with = cast_(get_r_type(x), with, x)); ++NP;

  R_xlen_t xn = vec_length(x);
  int where_size = vec_length(where);
  int with_size = vec_length(with);


#define CHECK_WHERE                                            \
  if (xi == 0) continue;                                       \
  if (!between<int>(xi, 1, xn)){                               \
    YIELD(NP);                                                 \
    Rf_error("`where` must be an integer vector of values between 1 and `length(x)`");\
  }

  R_xlen_t xi;
  R_xlen_t withi = 0;

  // Shallow copy lists and deep copy data of vectors
  if (!internal_in_place){
    if (Rf_isVectorList(x)){
      SHIELD(x = cpp_shallow_copy(x)); ++NP;
    } else {
      SHIELD(x = cpp_semi_copy(x)); ++NP;
    }
  }

  switch (get_r_type(x)){

  case r_null: {
    break;
  }
  case r_lgl:
  case r_int:
  case r_fct: {

    int* RESTRICT p_x = INTEGER(x);
    const int* RESTRICT p_with = INTEGER_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      p_x[xi - 1] = p_with[withi];
    }
    break;
  }
  case r_dbl: {

    double* RESTRICT p_x = REAL(x);
    const double* RESTRICT p_with = REAL_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      p_x[xi - 1] = p_with[withi];
    }
    break;
  }

  case r_int64: {

    int64_t* RESTRICT p_x = INTEGER64_PTR(x);
    const int64_t* RESTRICT p_with = INTEGER64_PTR_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      p_x[xi - 1] = p_with[withi];
    }
    break;
  }

  case r_chr: {

    const SEXP *p_with = STRING_PTR_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      SET_STRING_ELT(x, xi - 1, p_with[withi]);
    }
    break;
  }

  case r_cplx: {

    Rcomplex* p_x = COMPLEX(x);
    const Rcomplex *p_with = COMPLEX_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      --xi;
      p_x[xi].r = p_with[withi].r;
      p_x[xi].i = p_with[withi].i;
    }
    break;
  }

  case r_raw: {

    const Rbyte *p_with = RAW_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      SET_RAW_ELT(x, xi - 1, p_with[withi]);
    }
    break;
  }

  case r_date:
  case r_pxct: {
    SEXP x_cls = SHIELD(Rf_getAttrib(x, R_ClassSymbol)); ++NP;
    Rf_classgets(x, R_NilValue);
    static_cast<void>(cpp_assign(x, where, with, true));
    Rf_classgets(x, x_cls);
    break;
  }

  case r_list: {
    const SEXP *p_with = VECTOR_PTR_RO(with);

    for (int i = 0; i < where_size; withi = (++withi == with_size) ? 0 : withi, ++i){
      xi = p_where[i];
      CHECK_WHERE
      SET_VECTOR_ELT(x, xi - 1, p_with[withi]);
    }
    break;
  }

  case r_rcrd: {

    // Treat vctrs_rcrd as df and use r_df method

    SEXP x_cls = SHIELD(Rf_getAttrib(x, R_ClassSymbol)); ++NP;
    SEXP x_nms = SHIELD(get_names(x)); ++NP;
    SHIELD(x = cpp_list_as_df(x)); ++NP;
    SHIELD(with = cpp_list_as_df(with)); ++NP;
    SHIELD(x = cpp_assign(x, where, with, in_place)); ++NP;
    clear_attributes(x);
    set_names(x, x_nms);
    Rf_classgets(x, x_cls);
    break;
  }

  case r_df: {

    const SEXP *p_x = VECTOR_PTR_RO(x);
    const SEXP *p_with = VECTOR_PTR_RO(with);

    int ncol = Rf_length(x);

    if (ncol != Rf_length(with)){
      YIELD(NP);
      Rf_error("`ncol(x)` must equal `ncol(with)`");
    }

    SEXP x_names = SHIELD(get_names(x)); ++NP;
    SEXP with_names = SHIELD(get_names(with)); ++NP;

    if (!R_compute_identical(x_names, with_names, 0)){
      YIELD(NP);
      Rf_error("Column names must be identical between `x` and `with`");
    }

    for (int j = 0; j < ncol; ++j){
      SET_VECTOR_ELT(x, j, cpp_assign(p_x[j], where, p_with[j], in_place));
    }

    break;
  }

  case r_unk: {
    SHIELD(x = base_assign(x, where, with)); ++NP;
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
