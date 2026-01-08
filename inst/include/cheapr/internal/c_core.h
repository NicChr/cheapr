#ifndef CHEAPR_C_CORE_H
#define CHEAPR_C_CORE_H

// cheapr Core definitions and templates
// License: MIT

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_symbols.h>
#include <cheapr/internal/r_env.h>
#include <cheapr/internal/r_limits.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_nas.h>
#include <cheapr/internal/r_methods.h>
#include <cheapr/internal/r_attrs.h>
#include <cheapr/internal/r_rtype_coerce.h>
#include <cheapr/internal/r_vector.h>
#include <cheapr/internal/r_factor.h>
#include <cheapr/internal/r_list.h>
#include <cheapr/internal/r_coerce.h>
#include <cheapr/internal/r_make_vec.h>
#include <cheapr/internal/r_exprs.h>
#include <cheapr/internal/r_fns.h>
#include <cheapr/internal/r_math.h>
#include <optional>
#include <type_traits>


namespace cheapr {


// Functions

namespace altrep {
inline bool is_altrep(SEXP x){
  return ALTREP(x);
}
}

namespace df {

inline bool is_df(SEXP x){
  return attr::inherits1(x, "data.frame");
}

inline int nrow(SEXP x){
  return Rf_length(attr::get_attr(x, symbol::row_names_sym));
}
inline int ncol(SEXP x){
  return Rf_length(x);
}
inline SEXP new_row_names(int n){
  if (n > 0){
    auto out = SHIELD(vec::new_vector<r_int>(2));
    out.set(0, na::integer);
    out.set(1, -n);
    YIELD(1);
    return out;
  } else {
    return vec::new_vector<r_int>(0);
  }
}
inline void set_row_names(SEXP x, int n){
  SEXP row_names = SHIELD(new_row_names(n));
  attr::set_attr(x, symbol::row_names_sym, row_names);
  YIELD(1);
}
}

namespace vec {

inline r_size_t old_length(SEXP x){
  return Rf_xlength(x);
}

inline SEXP shallow_copy(SEXP x){
  return Rf_shallow_duplicate(x);
}

// Compact seq generator as ALTREP, same as `seq_len()`
inline SEXP compact_seq_len(r_size_t n){
  if (n < 0){
    Rf_error("`n` must be >= 0");
  }
  if (n == 0){
    return vec::new_vector<r_int>(0);
  }
  SEXP colon_fn = SHIELD(fn::find_pkg_fun(":", "base", false));
  SEXP out = SHIELD(fn::eval_fn(colon_fn, env::base_env, 1, n));
  YIELD(2);
  return out;
}

// r_lgl not bool because bool can't be NA
inline r_lgl all_whole_numbers(SEXP x, r_dbl tol_, bool na_rm_){

  r_size_t n = Rf_xlength(x);

  // Use r_lgl instead of bool as r_lgl can hold NA
  r_lgl out = r_true;
  r_size_t na_count = 0;

  switch ( internal::CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP:
  case internal::CHEAPR_INT64SXP: {
    break;
  }
  case REALSXP: {
    auto xvec = r_vec<r_dbl>(x);
    for (r_size_t i = 0; i < n; ++i) {
      out = static_cast<r_lgl>(math::is_whole_number(xvec.get(i), tol_));
      na_count += is_r_na(out);
      if (out == r_false){
        break;
      }
    }
    if (out == r_true && !na_rm_ && na_count > 0){
      out = na::logical;
    } else if (na_rm_ && na_count == n){
      out = r_true;
    }
    break;
  }
  default: {
    out = r_false;
    break;
  }
  }
  return out;
}
}

namespace internal {

inline void add_attrs(SEXP x, r_vec<r_sexp> attrs) {

  if (is_null(x)){
    Rf_error("Cannot add attributes to `NULL`");
  }

  if (is_null(attrs)){
    return;
  }

  int32_t NP = 0;

  auto names = SHIELD(r_vec<r_str>(attr::get_old_names(attrs))); ++NP;

  if (is_null(names)){
    YIELD(NP);
    Rf_error("attributes must be a named list");
  }

  r_sym attr_nm;

  int n = names.length();

  for (int i = 0; i < n; ++i){
    if (!(names.get(i) == blank_r_string)){
      attr_nm = as<r_sym>(names.get(i));
      if (address(x) == address(attrs.get(i))){
        SEXP dup_attr = SHIELD(Rf_duplicate(attrs.get(i))); ++NP;
        attr::set_attr(x, attr_nm, dup_attr);
      } else {
        attr::set_attr(x, attr_nm, attrs.get(i));
      }
    }
  }

YIELD(NP);
}

}

namespace attr {

inline void clear_attrs(SEXP x){
  
  auto attrs = SHIELD(static_cast<r_vec<r_sexp>>(get_attrs(x)));

  if (is_null(attrs)){
    YIELD(1);
    return;
  }
  auto names = SHIELD(static_cast<r_vec<r_str>>(attr::get_old_names(attrs)));

  int n = attrs.length();
  
  for (r_size_t i = 0; i < n; ++i){
    r_sym target_sym = as<r_sym>(names.get(i));
    set_attr(x, target_sym, r_null);
  }
  YIELD(2);
}

template<typename... Args>
inline void modify_attrs(SEXP x, Args... args) {
  auto attrs = SHIELD(make_list(args...));
  internal::add_attrs(x, attrs);
  YIELD(1);
}

inline void set_attrs(SEXP x, SEXP attrs){
  if (!is_null(x)){
    clear_attrs(x);
    internal::add_attrs(x, static_cast<r_vec<r_sexp>>(attrs));
  }
}

}

namespace internal {

// A cleaner lambda-based alternative to
// using the canonical switch(TYPEOF(x))
//
// Pass both the SEXP and an auto variable inside a lambda
// and visit_vector() will assign the auto variable to the correct vector
// Then simply deduce its type (via decltype) for further manipulation
// To be used in a lambda
// E.g. visit_vector(x, [&](auto x_vec) {})

// One must account for objects like `NULL` and non-vectors outwith this method

template <class F>
decltype(auto) visit_vector(SEXP x, F&& f) {
  switch (CHEAPR_TYPEOF(x)) {
  case LGLSXP:          return f(r_vec<r_lgl>(x));
  case INTSXP:          return f(r_vec<r_int>(x));
  case CHEAPR_INT64SXP: return f(r_vec<r_int64>(x));
  case REALSXP:         return f(r_vec<r_dbl>(x));
  case STRSXP:          return f(r_vec<r_str>(x));
  case VECSXP:          return f(r_vec<r_sexp>(x));
  case CPLXSXP:         return f(r_vec<r_cplx>(x));
  case RAWSXP:          return f(r_vec<r_raw>(x));
  default:              Rf_error("`x` must be a vector");
  }
}

}

namespace vec {
inline SEXP deep_copy(SEXP x){
  int32_t NP = 0;
  SEXP out = r_null;
  r_size_t n = Rf_xlength(x);

  if (!is_null(x)){

    if (Rf_isVector(x)){
      out = internal::visit_vector(x, [&](auto vec) -> SEXP {

        using r_t = decltype(vec);
        
        auto local_out = SHIELD(r_t(n)); ++NP;

        if constexpr (is<r_t, r_vec<r_sexp>>){
          for (r_size_t i = 0; i < n; ++i){
            local_out.set(i, deep_copy(vec.get(i)));
          }
        } else {
          r_copy_n(local_out, vec, 0, n);
        }

        return local_out;
      });
    } else {
      out = SHIELD(Rf_duplicate(x)); ++NP;
    }

    auto attrs = SHIELD(static_cast<r_vec<r_sexp>>(attr::get_attrs(x))); ++NP;
    int n_attrs = attrs.length();
    for (r_size_t i = 0; i < n_attrs; ++i){
      attrs.set(i, deep_copy(attrs.get(i)));
    }
    attr::set_attrs(out, attrs);
  }

  YIELD(NP);
  return out;
}

}


// We call R fn`cheapr::set_threads` to make sure the R option is set
inline void set_threads(uint16_t n){
  uint16_t max_threads = OMP_MAX_THREADS;
  uint16_t threads = std::min(n, max_threads);
  SEXP cheapr_set_threads = SHIELD(fn::find_pkg_fun("set_threads", "cheapr", true));
  SEXP r_threads = SHIELD(vec::as_vector(as<r_int>(threads)));
  SHIELD(fn::eval_fn(cheapr_set_threads, R_BaseEnv, r_threads));
  YIELD(3);
}


namespace internal {

// Wrap any callable f, and return a new callable that:
//   - takes (auto&&... args)
//   - calls f(args...) inside cpp11::unwind_protect

// Like cpp11::safe but works also  for variadic fns
template <typename F>
auto r_safe_impl(F f) {
  return [f](auto&&... args)
    -> decltype(f(std::forward<decltype(args)>(args)...)) {

      using result_t = decltype(f(std::forward<decltype(args)>(args)...));

      if constexpr (std::is_void_v<result_t>) {
        cpp11::unwind_protect([&] {
          f(std::forward<decltype(args)>(args)...);
        });
        // no return; result_t is void
      } else {
        return cpp11::unwind_protect([&]() -> result_t {
          return f(std::forward<decltype(args)>(args)...);
        });
      }
    };
}
}

#define r_safe(F)                                                                      \
internal::r_safe_impl(                                                                 \
  [&](auto&&... args)                                                                  \
    -> decltype(F(std::forward<decltype(args)>(args)...)) {                            \
      return F(std::forward<decltype(args)>(args)...);                                 \
    }                                                                                  \
)

} // End of cheapr namespace

#endif
