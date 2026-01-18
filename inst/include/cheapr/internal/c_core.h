#ifndef CHEAPR_C_CORE_H
#define CHEAPR_C_CORE_H

// cheapr Core definitions and templates
// License: MIT

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_symbols.h>
#include <cheapr/internal/r_env.h>
#include <cheapr/internal/r_limits.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_nas.h>
#include <cheapr/internal/r_methods.h>
#include <cheapr/internal/r_rtype_coerce.h>
#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_attrs.h>
#include <cheapr/internal/r_factor.h>
#include <cheapr/internal/r_list.h>
#include <cheapr/internal/r_coerce.h>
#include <cheapr/internal/r_make_vec.h>
// #include <cheapr/internal/r_df.h>
#include <cheapr/internal/r_exprs.h>
#include <cheapr/internal/r_fns.h>
#include <cheapr/internal/r_math.h>
#include <optional>
#include <type_traits>


namespace cheapr {


// Functions

namespace altrep {
inline bool is_altrep(r_sexp x){
  return ALTREP(x);
}
}

namespace attr {
template<typename... Args>
inline void modify_attrs(r_sexp x, Args... args) {
  auto attrs = make_vec<r_sexp>(args...);
  internal::modify_attrs_impl(x, attrs);
}
}

namespace vec {

inline r_size_t old_length(r_sexp x){
  return Rf_xlength(x);
}

inline r_sexp shallow_copy(r_sexp x){
  return r_sexp(Rf_shallow_duplicate(x)); 
}

// Compact seq generator as ALTREP, same as `seq_len()`
inline r_vec<r_int> compact_seq_len(r_size_t n){
  if (n < 0){
    cpp11::stop("`n` must be >= 0");
  }
  if (n == 0){
    return r_vec<r_int>();
  }
  r_sexp colon_fn = fn::find_pkg_fun(":", "base", false);
  r_sexp out = fn::eval_fn(colon_fn, env::base_env, 1, n);
  return r_vec<r_int>(out);
}

// r_lgl not bool because bool can't be NA
inline r_lgl all_whole_numbers(r_sexp x, r_dbl tol_, bool na_rm_){

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
    auto xvec = as<r_vec<r_dbl>>(x);
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

namespace vec {
inline r_sexp deep_copy(r_sexp x){
  r_sexp out = r_null;
  r_size_t n = Rf_xlength(x);

  if (!x.is_null()){
    if (Rf_isVector(x)){
      out = internal::visit_vector(x, [&](auto vec) -> r_sexp {

        using r_t = decltype(vec);

        auto local_out = r_t(n);

        if constexpr (is<r_t, r_vec<r_sexp>>){
          for (r_size_t i = 0; i < n; ++i){
            local_out.set(i, deep_copy(vec.get(i)));
          }
        } else {
          r_copy_n(local_out, vec, 0, n);
        }

        return local_out.sexp;
      });
    } else {
      out = r_sexp(Rf_duplicate(x));
    }

    auto attrs = attr::get_attrs(x);
    int n_attrs = attrs.length();
    for (int i = 0; i < n_attrs; ++i){
      attrs.set(i, deep_copy(attrs.get(i)));
    }
    attr::set_attrs(out, attrs);
  }

  return out;
}

}


// We call R fn`cheapr::set_threads` to make sure the R option is set
inline void set_threads(uint16_t n){
  uint16_t max_threads = OMP_MAX_THREADS;
  uint16_t threads = std::min(n, max_threads);
  r_sexp cheapr_set_threads = fn::find_pkg_fun("set_threads", "cheapr", true);
  auto r_threads = as_vector(as<r_int>(threads));
  fn::eval_fn(cheapr_set_threads, env::base_env, r_threads);
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
