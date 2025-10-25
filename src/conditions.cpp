#include "cheapr.h"

// Fast SIMD vectorised if-else

SEXP if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){

  int yes_type = TYPEOF(yes);

  if (TYPEOF(condition) != LGLSXP){
    Rf_error("condition must be a logical vector");
  }
  if (yes_type != TYPEOF(no)){
    Rf_error("`typeof(yes)` must match `typeof(no)`");
  }
  if (yes_type != TYPEOF(na)){
    Rf_error("`typeof(yes)` must match `typeof(na)`");
  }
  R_xlen_t n = Rf_xlength(condition);
  R_xlen_t yes_size = Rf_xlength(yes);
  R_xlen_t no_size = Rf_xlength(no);
  R_xlen_t na_size = Rf_xlength(na);

  if (yes_size != 1 && yes_size != n){
    Rf_error("`length(yes)` must be 1 or `length(condition)`");
  }
  if (no_size != 1 && no_size != n){
    Rf_error("`length(no)` must be 1 or `length(condition)`");
  }
  if (na_size != 1 && na_size != n){
    Rf_error("`length(na)` must be 1 or `length(condition)`");
  }

  bool yes_scalar = yes_size == 1;
  bool no_scalar = no_size == 1;
  bool na_scalar = na_size == 1;
  bool all_scalar = yes_scalar && no_scalar && na_scalar;
  int lgl;

#define VECTORISED_IF_ELSE                                         \
  for (R_xlen_t i = 0; i < n; ++i){                                \
    switch(p_x[i]){                                                \
    case 1: {                                                      \
  p_out[i] = p_yes[yes_scalar ? 0 : i];                            \
  break;                                                           \
}                                                                  \
    case 0: {                                                      \
      p_out[i] = p_no[no_scalar ? 0 : i];                          \
      break;                                                       \
    }                                                              \
    default: {                                                     \
      p_out[i] = p_na[na_scalar ? 0 : i];                          \
      break;                                                       \
    }                                                              \
    }                                                              \
  }                                                                \

#define SCALAR_IF_ELSE                                               \
  for (R_xlen_t i = 0; i < n; ++i){                                  \
    lgl = p_x[i];                                                    \
    if (lgl == 1){                                                   \
      p_out[i] = yes_value;                                          \
    } else if (lgl == 0){                                            \
      p_out[i] = no_value;                                           \
    } else {                                                         \
      p_out[i] = na_value;                                           \
    }                                                                \
  }                                                                  \

  const int* RESTRICT p_x = INTEGER_RO(condition);
  SEXP out = SHIELD(new_vec(yes_type, n));

  switch (yes_type){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int* RESTRICT p_out = INTEGER(out);
    const int *p_yes = INTEGER(yes);
    const int *p_no = INTEGER(no);
    const int *p_na = INTEGER(na);

    if (all_scalar){
      const int yes_value = p_yes[0];
      const int no_value = p_no[0];
      const int na_value = p_na[0];
      OMP_FOR_SIMD
      SCALAR_IF_ELSE
    } else {
      VECTORISED_IF_ELSE
    }
    break;
  }
  case REALSXP: {
    double* RESTRICT p_out = REAL(out);
    const double *p_yes = REAL(yes);
    const double *p_no = REAL(no);
    const double *p_na = REAL(na);

    if (all_scalar){
      const double yes_value = p_yes[0];
      const double no_value = p_no[0];
      const double na_value = p_na[0];
      OMP_FOR_SIMD
      SCALAR_IF_ELSE
    } else {
      VECTORISED_IF_ELSE
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_yes = STRING_PTR_RO(yes);
    const SEXP *p_no = STRING_PTR_RO(no);
    const SEXP *p_na = STRING_PTR_RO(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case 1: {
      SET_STRING_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case 0: {
        SET_STRING_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_STRING_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  case CPLXSXP: {
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

          if (lgl == 1){
            p_out[i].r = yes_value_re;
            p_out[i].i = yes_value_im;
          } else if (lgl == 0){
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
        case 1: {
        SET_COMPLEX_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
        break;
      }
        case 0: {
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
  case RAWSXP: {
    const Rbyte *p_yes = RAW(yes);
    const Rbyte *p_no = RAW(no);
    const Rbyte *p_na = RAW(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case 1: {
      SET_RAW_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case 0: {
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
  case VECSXP: {
    const SEXP *p_yes = VECTOR_PTR_RO(yes);
    const SEXP *p_no = VECTOR_PTR_RO(no);
    const SEXP *p_na = VECTOR_PTR_RO(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case 1: {
      SET_VECTOR_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case 0: {
        SET_VECTOR_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_VECTOR_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  default: {
    YIELD(1);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(yes_type));
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){

  int32_t NP = 0;

  if (is_null(na)){
    SHIELD(na = cpp_na_init(no, 1)); ++NP;
  }

  SEXP args = SHIELD(new_vec(VECSXP, 3)); ++NP;
  SET_VECTOR_ELT(args, 0, yes);
  SET_VECTOR_ELT(args, 1, no);
  SET_VECTOR_ELT(args, 2, na);
  SEXP out = R_NilValue;

  // Fast method for bare atomic vectors
  if (is_bare_atomic(yes) && is_bare_atomic(no) && is_bare_atomic(na)){
    SHIELD(args = cpp_cast(args)); ++NP;
    const SEXP *p_args = VECTOR_PTR_RO(args);

    yes = p_args[0];
    no = p_args[1];
    na = p_args[2];

    SHIELD(out = if_else(condition, yes, no, na)); ++NP;
    YIELD(NP);
    return out;
  }

  // Below is a catch-all method

  SEXP templates = SHIELD(new_vec(VECSXP, 3)); ++NP;
  SET_VECTOR_ELT(templates, 0, slice_loc(yes, 1));
  SET_VECTOR_ELT(templates, 1, slice_loc(no, 1));
  SET_VECTOR_ELT(templates, 2, slice_loc(na, 1));

  if (!Rf_isLogical(condition)){
    YIELD(NP);
    Rf_error("condition must be a logical vector");
  }
  if (vec_length(yes) != 1 && vec_length(yes) != Rf_xlength(condition)){
    YIELD(NP);
    Rf_error("`length(true)` must be 1 or `length(condition)`");
  }
  if (vec_length(no) != 1 && vec_length(no) != Rf_xlength(condition)){
    YIELD(NP);
    Rf_error("`length(false)` must be 1 or `length(condition)`");
  }
  if (vec_length(na) != 1 && vec_length(na) != Rf_xlength(condition)){
    YIELD(NP);
    Rf_error("`length(na)` must be 1 or `length(condition)`");
  }

  SEXP assign_sym = SHIELD(install_utf8("[<-")); ++NP;
  SEXP cast_expr = SHIELD(
    Rf_lang3(
      R_TripleColonSymbol, install_utf8("cheapr"), install_utf8("cast")
    )
  ); ++NP;
  SEXP cast_fn = SHIELD(Rf_eval(cast_expr, R_BaseEnv)); ++NP;
  SEXP cast_template = R_NilValue;
  SEXP expr = R_NilValue;

  SHIELD(cast_template = cpp_c(templates)); ++NP;
  SHIELD(cast_template = slice_loc(cast_template, 0)); ++NP;
  SHIELD(expr = Rf_lang3(cast_fn, yes, cast_template)); ++NP;
  SHIELD(yes = Rf_eval(expr, R_GetCurrentEnv())); ++NP;
  SHIELD(expr = Rf_lang3(cast_fn, no, cast_template)); ++NP;
  SHIELD(no = Rf_eval(expr, R_GetCurrentEnv())); ++NP;
  SHIELD(expr = Rf_lang3(cast_fn, na, cast_template)); ++NP;
  SHIELD(na = Rf_eval(expr, R_GetCurrentEnv())); ++NP;

  SEXP lgl_val_counts = SHIELD(cpp_lgl_count(condition)); ++NP;
  SHIELD(lgl_val_counts = coerce_vec(lgl_val_counts, REALSXP)); ++NP;
  R_xlen_t n_true = REAL_RO(lgl_val_counts)[0];
  R_xlen_t n_false = REAL_RO(lgl_val_counts)[1];
  R_xlen_t n_na = REAL_RO(lgl_val_counts)[2];

  if (n_true == Rf_xlength(condition)){
    if (Rf_xlength(yes) == 1){
      SHIELD(out = cpp_rep_len(yes, Rf_xlength(condition))); ++NP;
      YIELD(NP);
      return out;
    } else {
      YIELD(NP);
      return yes;
    }
  }

  if (n_false == Rf_xlength(condition)){
    if (Rf_xlength(no) == 1){
      SHIELD(out = cpp_rep_len(no, Rf_xlength(condition))); ++NP;
      YIELD(NP);
      return out;
    } else {
      YIELD(NP);
      return no;
    }
  }

  if (n_na == Rf_xlength(condition)){
    if (Rf_xlength(na) == 1){
      SHIELD(out = cpp_rep_len(na, Rf_xlength(condition))); ++NP;
      YIELD(NP);
      return out;
    } else {
      YIELD(NP);
      return na;
    }
  }

  // Assuming the else part is most likely to be most prominent
  if (vec_length(no) == Rf_xlength(condition)){
    out = no;
  } else {
    SHIELD(out = cpp_rep_len(no, Rf_xlength(condition))); ++NP;
  }

  SEXP lgl_locs = SHIELD(cpp_lgl_locs(
    condition, n_true, n_false, true, false, true
  )); ++NP;
  SEXP true_locs = get_list_element(lgl_locs, make_utf8_char("true"));
  SEXP na_locs = get_list_element(lgl_locs, make_utf8_char("na"));


  SEXP replace = R_NilValue;

  if (vec_length(yes) == 1){
    replace = yes;
  } else {
    SHIELD(replace = cpp_sset(yes, true_locs, true)); ++NP;
  }

  SHIELD(expr = Rf_lang4(assign_sym, out, true_locs, replace)); ++NP;
  SEXP arg = CDR(expr);
  SET_TAG(arg, install_utf8("x"));
  arg = CDR(arg);
  SET_TAG(arg, install_utf8("i"));
  arg = CDR(arg);
  SET_TAG(arg, install_utf8("value"));

  SHIELD(out = Rf_eval(expr, R_GetCurrentEnv())); ++NP;

  if (vec_length(yes) == 1){
    replace = na;
  } else {
    SHIELD(replace = cpp_sset(na, na_locs, true)); ++NP;
  }
  SHIELD(expr = Rf_lang4(assign_sym, out, na_locs, replace)); ++NP;
  arg = CDR(expr);
  SET_TAG(arg, install_utf8("x"));
  arg = CDR(arg);
  SET_TAG(arg, install_utf8("i"));
  arg = CDR(arg);
  SET_TAG(arg, install_utf8("value"));
  SHIELD(out = Rf_eval(expr, R_GetCurrentEnv())); ++NP;
  YIELD(NP);
  return out;
}


