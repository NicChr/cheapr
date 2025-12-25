#include "cheapr.h"

#define CHEAPR_VECTORISED_IF_ELSE(OMP_ROUTINE)                   \
if (all_not_scalar){                                             \
  OMP_ROUTINE                                                    \
  for (R_xlen_t i = 0; i < n; ++i){                              \
    if (p_x[i] == r_true){                                       \
      set_value(p_out, i, p_yes[i]);                               \
    } else if (p_x[i] == r_false){                               \
      set_value(p_out, i, p_no[i]);                                \
    } else {                                                     \
      set_value(p_out, i, p_na[i]);                                \
    }                                                            \
  }                                                              \
} else {                                                         \
  OMP_ROUTINE                                                    \
  for (R_xlen_t i = 0; i < n; ++i){                              \
    if (p_x[i] == r_true){                                       \
      set_value(p_out, i, p_yes[yes_scalar ? 0 : i]);              \
    } else if (p_x[i] == r_false){                               \
      set_value(p_out, i, p_no[no_scalar ? 0 : i]);                \
    } else {                                                     \
      set_value(p_out, i, p_na[na_scalar ? 0 : i]);                \
    }                                                            \
  }                                                              \
}


#define CHEAPR_SCALAR_IF_ELSE                                                      \
for (R_xlen_t i = 0; i < n; ++i){                                                  \
  lgl = p_x[i];                                                                    \
  if (lgl == r_true){                                                              \
    set_value(p_out, i, yes_value);                                                  \
  } else if (lgl == r_false){                                                      \
    set_value(p_out, i, no_value);                                                   \
  } else {                                                                         \
    set_value(p_out, i, na_value);                                                   \
  }                                                                                \
}

#define CHEAPR_DO_IF_ELSE                                                 \
if (all_scalar){                                                          \
  if (n_threads > 1){                                                     \
    OMP_PARALLEL_FOR_SIMD(n_threads)                                      \
    CHEAPR_SCALAR_IF_ELSE                                                 \
  } else {                                                                \
    OMP_SIMD                                                              \
    CHEAPR_SCALAR_IF_ELSE                                                 \
  }                                                                       \
} else {                                                                  \
  if (n_threads > 1){                                                     \
    CHEAPR_VECTORISED_IF_ELSE(OMP_PARALLEL_FOR_SIMD(n_threads))           \
  } else {                                                                \
    CHEAPR_VECTORISED_IF_ELSE(OMP_SIMD)                                   \
  }                                                                       \
}



// Fast SIMD vectorised if-else

[[cpp11::register]]
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){

  if (!Rf_isLogical(condition)){
    Rf_error("condition must be a logical vector");
  }

  const r_bool_t* RESTRICT p_x = logical_ptr_ro(condition);

  int32_t NP = 0;

  if (is_null(na)){
    SHIELD(na = cpp_na_init(no, 1)); ++NP;
  }

  SEXP args = SHIELD(make_list(yes, no, na)); ++NP;
  SHIELD(args = cpp_cast_common(args)); ++NP;
  yes = VECTOR_ELT(args, 0);
  no = VECTOR_ELT(args, 1);
  na = VECTOR_ELT(args, 2);

  R_xlen_t n = vec::length(condition);

  if (n == 0){
    YIELD(NP);
    return yes;
  }

  r_type common = get_r_type(yes);

  if (common != R_unk){

    SEXP out = r_null;

    R_xlen_t yes_size = vec::length(yes);
    R_xlen_t no_size = vec::length(no);
    R_xlen_t na_size = vec::length(na);

    if (yes_size != 1 && yes_size != n){
      YIELD(NP);
      Rf_error("`vector_length(true)` must be 1 or `length(condition)`");
    }
    if (no_size != 1 && no_size != n){
      YIELD(NP);
      Rf_error("`vector_length(false)` must be 1 or `length(condition)`");
    }
    if (na_size != 1 && na_size != n){
      YIELD(NP);
      Rf_error("`vector_length(na)` must be 1 or `length(condition)`");
    }

    bool yes_scalar = yes_size == 1;
    bool no_scalar = no_size == 1;
    bool na_scalar = na_size == 1;
    bool all_scalar = yes_scalar && no_scalar && na_scalar;
    bool all_not_scalar = !yes_scalar && !no_scalar && !na_scalar;
    r_bool_t lgl;

    int n_threads = calc_threads(n);

    switch (common){
    case R_null: {
      break;
    }
    case R_lgl: {
      SHIELD(out = init<r_logicals_t>(n, false)); ++NP;
      r_bool_t* RESTRICT p_out = logical_ptr(out);
      const r_bool_t *p_yes = logical_ptr_ro(yes);
      const r_bool_t *p_no = logical_ptr_ro(no);
      const r_bool_t *p_na = logical_ptr_ro(na);

      const r_bool_t yes_value = p_yes[0];
      const r_bool_t no_value = p_no[0];
      const r_bool_t na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_int: {
      SHIELD(out = init<r_integers_t>(n, false)); ++NP;
      int* RESTRICT p_out = integer_ptr(out);
      const int *p_yes = integer_ptr_ro(yes);
      const int *p_no = integer_ptr_ro(no);
      const int *p_na = integer_ptr_ro(na);
      const int yes_value = p_yes[0];
      const int no_value = p_no[0];
      const int na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_int64: {
      SHIELD(out = init<r_integers64_t>(n, false)); ++NP;
      int64_t* RESTRICT p_out = INTEGER64_PTR(out);
      const int64_t *p_yes = integer64_ptr_ro(yes);
      const int64_t *p_no = integer64_ptr_ro(no);
      const int64_t *p_na = integer64_ptr_ro(na);
      const int64_t yes_value = p_yes[0];
      const int64_t no_value = p_no[0];
      const int64_t na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_dbl: {
      SHIELD(out = init<r_doubles_t>(n, false)); ++NP;
      double* RESTRICT p_out = real_ptr(out);
      const double *p_yes = real_ptr_ro(yes);
      const double *p_no = real_ptr_ro(no);
      const double *p_na = real_ptr_ro(na);
      const double yes_value = p_yes[0];
      const double no_value = p_no[0];
      const double na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_chr: {

      SHIELD(out = init<r_characters_t>(n, false)); ++NP;
      SEXP p_out = out;

      const r_string_t *p_yes = string_ptr_ro(yes);
      const r_string_t *p_no = string_ptr_ro(no);
      const r_string_t *p_na = string_ptr_ro(na);

      if (all_scalar){
        const r_string_t yes_value = p_yes[0];
        const r_string_t no_value = p_no[0];
        const r_string_t na_value = p_na[0];
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE(OMP_DO_NOTHING)
      }
      break;
    }
    case R_cplx: {

      SHIELD(out = init<r_complexes_t>(n, false)); ++NP;
      r_complex_t *p_out = complex_ptr(out);

      const r_complex_t *p_yes = complex_ptr(yes);
      const r_complex_t *p_no = complex_ptr(no);
      const r_complex_t *p_na = complex_ptr(na);
      const r_complex_t yes_value = p_yes[0];
      const r_complex_t no_value = p_no[0];
      const r_complex_t na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_raw: {

      SHIELD(out = init<r_raws_t>(n, false)); ++NP;
      r_byte_t *p_out = raw_ptr(out);

      const r_byte_t *p_yes = raw_ptr(yes);
      const r_byte_t *p_no = raw_ptr(no);
      const r_byte_t *p_na = raw_ptr(na);
      const r_byte_t yes_value = p_yes[0];
      const r_byte_t no_value = p_no[0];
      const r_byte_t na_value = p_na[0];
      CHEAPR_DO_IF_ELSE

      break;
    }
    case R_fct: {
      SHIELD(out = cpp_na_init(yes, n)); ++NP;
      int* RESTRICT p_out = integer_ptr(out);

      const int *p_yes = integer_ptr_ro(yes);
      const int *p_no = integer_ptr_ro(no);
      const int *p_na = integer_ptr_ro(na);
      const int yes_value = p_yes[0];
      const int no_value = p_no[0];
      const int na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_date: {

      SHIELD(out = cpp_na_init(yes, n)); ++NP;

      if (TYPEOF(out) == INTSXP){
        int* RESTRICT p_out = integer_ptr(out);
        const int *p_yes = integer_ptr_ro(yes);
        const int *p_no = integer_ptr_ro(no);
        const int *p_na = integer_ptr_ro(na);
        const int yes_value = p_yes[0];
        const int no_value = p_no[0];
        const int na_value = p_na[0];
        CHEAPR_DO_IF_ELSE
      } else {
        double* RESTRICT p_out = real_ptr(out);
        const double *p_yes = real_ptr_ro(yes);
        const double *p_no = real_ptr_ro(no);
        const double *p_na = real_ptr_ro(na);
        const double yes_value = p_yes[0];
        const double no_value = p_no[0];
        const double na_value = p_na[0];
        CHEAPR_DO_IF_ELSE
      }

      break;
    }
    case R_pxt: {

      SHIELD(out = cpp_na_init(yes, n)); ++NP;
      double* RESTRICT p_out = real_ptr(out);

      const double *p_yes = real_ptr_ro(yes);
      const double *p_no = real_ptr_ro(no);
      const double *p_na = real_ptr_ro(na);
      const double yes_value = p_yes[0];
      const double no_value = p_no[0];
      const double na_value = p_na[0];
      CHEAPR_DO_IF_ELSE
      break;
    }
    case R_list: {

      SHIELD(out = init<r_list_t>(n, false)); ++NP;
      SEXP p_out = out;

      const SEXP *p_yes = list_ptr_ro(yes);
      const SEXP *p_no = list_ptr_ro(no);
      const SEXP *p_na = list_ptr_ro(na);

      if (all_scalar){

        const SEXP yes_value = p_yes[0];
        const SEXP no_value = p_no[0];
        const SEXP na_value = p_na[0];
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE(OMP_DO_NOTHING)
      }
      break;
    }
    case R_df: {

      SHIELD(out = cpp_na_init(yes, n)); ++NP;

      int ncol = df::ncol(yes);

      if (ncol != df::ncol(no)){
        YIELD(NP);
        Rf_error("`ncol(yes)` must equal `ncol(no)`");
      }
      if (ncol != df::ncol(na)){
        YIELD(NP);
        Rf_error("`ncol(yes)` must equal `ncol(na)`");
      }

      SEXP col_names = SHIELD(attr::get_old_names(yes)); ++NP;
      SEXP col_names2 = SHIELD(attr::get_old_names(no)); ++NP;
      SEXP col_names3 = SHIELD(attr::get_old_names(na)); ++NP;

      if (!R_compute_identical(col_names, col_names2, 0)){
        YIELD(NP);
        Rf_error("Column names must be identical between `yes` and `no`");
      }
      if (!R_compute_identical(col_names, col_names3, 0)){
        YIELD(NP);
        Rf_error("Column names must be identical between `yes` and `na`");
      }

      for (int j = 0; j < ncol; ++j){
        SET_VECTOR_ELT(out, j, cpp_if_else(condition, VECTOR_ELT(yes, j), VECTOR_ELT(no, j), VECTOR_ELT(na, j)));
      }
      break;
    }
    default: {
      YIELD(NP);
      Rf_error("%s cannot handle an object of type %s", __func__, r_type_char(yes));
    }
    }
    YIELD(NP);
    return out;
  } else {
    // We're calling an R function instead of doing it here in C/C++
    // to take advantage of the fact that `[<-` avoids creating copies due to
    // correct reference-tracking in R
    // If we call `[<-` directly then unnecessary copies are made
    SEXP out = SHIELD(eval_pkg_fun("if_else2", "cheapr", env::base_env, yes, no, na)); ++NP;
    YIELD(NP);
    return out;
  }
}


