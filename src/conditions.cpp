#include "cheapr.h"

#define CHEAPR_VECTORISED_IF_ELSE                                    \
for (R_xlen_t i = 0; i < n; ++i){                                    \
  switch(p_x[i]){                                                    \
  case r_true: {                                                     \
    p_out[i] = p_yes[yes_scalar ? 0 : i];                            \
    break;                                                           \
  }                                                                  \
  case r_false: {                                                    \
    p_out[i] = p_no[no_scalar ? 0 : i];                              \
    break;                                                           \
  }                                                                  \
  default: {                                                         \
    p_out[i] = p_na[na_scalar ? 0 : i];                              \
    break;                                                           \
  }                                                                  \
  }                                                                  \
}

#define CHEAPR_SCALAR_IF_ELSE                                              \
for (R_xlen_t i = 0; i < n; ++i){                                          \
  lgl = p_x[i];                                                            \
  if (lgl == r_true){                                                      \
    p_out[i] = yes_value;                                                  \
  } else if (lgl == r_false){                                              \
    p_out[i] = no_value;                                                   \
  } else {                                                                 \
    p_out[i] = na_value;                                                   \
  }                                                                        \
}

// Fast SIMD vectorised if-else

[[cpp11::register]]
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){

  if (!Rf_isLogical(condition)){
    Rf_error("condition must be a logical vector");
  }

  int32_t NP = 0;

  if (is_null(na)){
    SHIELD(na = cpp_na_init(no, 1)); ++NP;
  }

  SEXP args = SHIELD(new_r_list(yes, no, na)); ++NP;
  SHIELD(args = cpp_cast_common(args)); ++NP;
  yes = VECTOR_ELT(args, 0);
  no = VECTOR_ELT(args, 1);
  na = VECTOR_ELT(args, 2);
  r_type common = get_r_type(yes);

  if (common != r_unk){
    SEXP out = R_NilValue;

    R_xlen_t n = vector_length(condition);
    R_xlen_t yes_size = vector_length(yes);
    R_xlen_t no_size = vector_length(no);
    R_xlen_t na_size = vector_length(na);

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
    r_boolean lgl;

    const r_boolean* RESTRICT p_x = BOOLEAN_RO(condition);

    switch (common){
    case r_null: {
      break;
    }
    case r_lgl: {
      SHIELD(out = init<r_logical_t>(n, false)); ++NP;
      r_boolean* RESTRICT p_out = BOOLEAN(out);
      const r_boolean *p_yes = BOOLEAN_RO(yes);
      const r_boolean *p_no = BOOLEAN_RO(no);
      const r_boolean *p_na = BOOLEAN_RO(na);

      if (all_scalar){
        const r_boolean yes_value = p_yes[0];
        const r_boolean no_value = p_no[0];
        const r_boolean na_value = p_na[0];
        OMP_FOR_SIMD
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE
      }
      break;
    }
    case r_int: {
      SHIELD(out = init<r_integer_t>(n, false)); ++NP;
      int* RESTRICT p_out = INTEGER(out);
      const int *p_yes = INTEGER_RO(yes);
      const int *p_no = INTEGER_RO(no);
      const int *p_na = INTEGER_RO(na);

      if (all_scalar){
        const int yes_value = p_yes[0];
        const int no_value = p_no[0];
        const int na_value = p_na[0];
        OMP_FOR_SIMD
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE
      }
      break;
    }
    case r_int64: {
      SHIELD(out = init<r_integer64_t>(n, false)); ++NP;
      int64_t* RESTRICT p_out = INTEGER64_PTR(out);
      const int64_t *p_yes = INTEGER64_PTR_RO(yes);
      const int64_t *p_no = INTEGER64_PTR_RO(no);
      const int64_t *p_na = INTEGER64_PTR_RO(na);

      if (all_scalar){
        const int64_t yes_value = p_yes[0];
        const int64_t no_value = p_no[0];
        const int64_t na_value = p_na[0];
        OMP_FOR_SIMD
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE
      }
      break;
    }
    case r_dbl: {
      SHIELD(out = init<r_numeric_t>(n, false)); ++NP;
      double* RESTRICT p_out = REAL(out);
      const double *p_yes = REAL_RO(yes);
      const double *p_no = REAL_RO(no);
      const double *p_na = REAL_RO(na);

      if (all_scalar){
        const double yes_value = p_yes[0];
        const double no_value = p_no[0];
        const double na_value = p_na[0];
        OMP_FOR_SIMD
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE
      }
      break;
    }
    case r_chr: {

      SHIELD(out = init<r_character_t>(n, false)); ++NP;

      const SEXP *p_yes = STRING_PTR_RO(yes);
      const SEXP *p_no = STRING_PTR_RO(no);
      const SEXP *p_na = STRING_PTR_RO(na);

      if (all_scalar){

        const SEXP yes_value = p_yes[0];
        const SEXP no_value = p_no[0];
        const SEXP na_value = p_na[0];

        for (R_xlen_t i = 0; i < n; ++i){
          lgl = p_x[i];

          if (lgl == r_true){
            SET_STRING_ELT(out, i, yes_value);
          } else if (lgl == r_false){
            SET_STRING_ELT(out, i, no_value);
          } else {
            SET_STRING_ELT(out, i, na_value);
          }
        }
      } else {

        for (R_xlen_t i = 0; i < n; ++i){
          switch(p_x[i]){
          case r_true: {
            SET_STRING_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
            break;
          }
          case r_false: {
            SET_STRING_ELT(out, i, p_no[no_scalar ? 0 : i]);
            break;
          }
          default: {
            SET_STRING_ELT(out, i, p_na[na_scalar ? 0 : i]);
            break;
          }
          }
        }
      }
      break;
    }
    case r_cplx: {

      SHIELD(out = init<r_complex_t>(n, false)); ++NP;
      Rcomplex *p_out = COMPLEX(out);

      const Rcomplex *p_yes = COMPLEX(yes);
      const Rcomplex *p_no = COMPLEX(no);
      const Rcomplex *p_na = COMPLEX(na);

      if (all_scalar){

        const double yes_value_re = p_yes[0].r;
        const double yes_value_im = p_yes[0].i;
        const double no_value_re = p_no[0].r;
        const double no_value_im = p_no[0].i;
        const double na_value_re = p_na[0].r;
        const double na_value_im = p_na[0].i;

        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n; ++i){
          lgl = p_x[i];

          if (lgl == r_true){
            p_out[i].r = yes_value_re;
            p_out[i].i = yes_value_im;
          } else if (lgl == r_false){
            p_out[i].r = no_value_re;
            p_out[i].i = no_value_im;
          } else {
            p_out[i].r = na_value_re;
            p_out[i].i = na_value_im;
          }
        }
      } else {

        for (R_xlen_t i = 0; i < n; ++i){
          switch(p_x[i]){
          case r_true: {
        SET_COMPLEX_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
        break;
      }
          case r_false: {
            SET_COMPLEX_ELT(out, i, p_no[no_scalar ? 0 : i]);
            break;
          }
          default: {
            SET_COMPLEX_ELT(out, i, p_na[na_scalar ? 0 : i]);
            break;
          }
          }
        }
      }
      break;
    }
    case r_raw: {

      SHIELD(out = init<r_raw_t>(n, false)); ++NP;

      const Rbyte *p_yes = RAW(yes);
      const Rbyte *p_no = RAW(no);
      const Rbyte *p_na = RAW(na);

      for (R_xlen_t i = 0; i < n; ++i){
        switch(p_x[i]){
        case r_true: {
        SET_RAW_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
        break;
      }
        case r_false: {
          SET_RAW_ELT(out, i, p_no[no_scalar ? 0 : i]);
          break;
        }
        default: {
          SET_RAW_ELT(out, i, p_na[na_scalar ? 0 : i]);
          break;
        }
        }
      }
      break;
    }
    case r_fct: {
      SHIELD(out = cpp_na_init(yes, n)); ++NP;
      int* RESTRICT p_out = INTEGER(out);

      const int *p_yes = INTEGER_RO(yes);
      const int *p_no = INTEGER_RO(no);
      const int *p_na = INTEGER_RO(na);

      if (all_scalar){
        const int yes_value = p_yes[0];
        const int no_value = p_no[0];
        const int na_value = p_na[0];
        OMP_FOR_SIMD
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE
      }
      break;
    }
    case r_date: {

      SHIELD(out = cpp_na_init(yes, n)); ++NP;

      if (TYPEOF(out) == INTSXP){
        int* RESTRICT p_out = INTEGER(out);
        const int *p_yes = INTEGER_RO(yes);
        const int *p_no = INTEGER_RO(no);
        const int *p_na = INTEGER_RO(na);

        if (all_scalar){
          const int yes_value = p_yes[0];
          const int no_value = p_no[0];
          const int na_value = p_na[0];
          OMP_FOR_SIMD
          CHEAPR_SCALAR_IF_ELSE
        } else {
          CHEAPR_VECTORISED_IF_ELSE
        }
      } else {
        double* RESTRICT p_out = REAL(out);
        const double *p_yes = REAL_RO(yes);
        const double *p_no = REAL_RO(no);
        const double *p_na = REAL_RO(na);

        if (all_scalar){
          const double yes_value = p_yes[0];
          const double no_value = p_no[0];
          const double na_value = p_na[0];
          OMP_FOR_SIMD
          CHEAPR_SCALAR_IF_ELSE
        } else {
          CHEAPR_VECTORISED_IF_ELSE
        }
      }

      break;
    }
    case r_pxct: {

      SHIELD(out = cpp_na_init(yes, n)); ++NP;
      double* RESTRICT p_out = REAL(out);

      const double *p_yes = REAL_RO(yes);
      const double *p_no = REAL_RO(no);
      const double *p_na = REAL_RO(na);

      if (all_scalar){
        const double yes_value = p_yes[0];
        const double no_value = p_no[0];
        const double na_value = p_na[0];
        OMP_FOR_SIMD
        CHEAPR_SCALAR_IF_ELSE
      } else {
        CHEAPR_VECTORISED_IF_ELSE
      }
      break;
    }
    case r_list: {

      SHIELD(out = init<r_list_t>(n, false)); ++NP;

      const SEXP *p_yes = LIST_PTR_RO(yes);
      const SEXP *p_no = LIST_PTR_RO(no);
      const SEXP *p_na = LIST_PTR_RO(na);

      if (all_scalar){

        const SEXP yes_value = p_yes[0];
        const SEXP no_value = p_no[0];
        const SEXP na_value = p_na[0];

        for (R_xlen_t i = 0; i < n; ++i){
          lgl = p_x[i];

          if (lgl == r_true){
            SET_VECTOR_ELT(out, i, yes_value);
          } else if (lgl == r_false){
            SET_VECTOR_ELT(out, i, no_value);
          } else {
            SET_VECTOR_ELT(out, i, na_value);
          }
        }
      } else {

        for (R_xlen_t i = 0; i < n; ++i){
          switch(p_x[i]){
          case r_true: {
            SET_VECTOR_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
            break;
          }
          case r_false: {
            SET_VECTOR_ELT(out, i, p_no[no_scalar ? 0 : i]);
            break;
          }
          default: {
            SET_VECTOR_ELT(out, i, p_na[na_scalar ? 0 : i]);
            break;
          }
          }
        }
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
    SEXP if_else_fn = SHIELD(find_pkg_fun("if_else2", "cheapr", true)); ++NP;

    // We're calling an R function instead of doing it here in C/C++
    // to take advantage of the fact that `[<-` avoids creating copies due to
    // correct reference-tracking in R
    // If we call `[<-` directly then unnecessary copies are made

    SEXP expr = SHIELD(Rf_lang5(if_else_fn, condition, yes, no, na)); ++NP;
    SEXP out = SHIELD(Rf_eval(expr, R_GetCurrentEnv())); ++NP;
    YIELD(NP);
    return out;
  }
}


