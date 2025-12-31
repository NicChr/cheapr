#ifndef CHEAPR_R_VECTOR
#define CHEAPR_R_VECTOR

#include <cpp11.hpp>
#include <r_types.h>
#include <r_concepts.h>

namespace cheapr {

namespace internal {

inline SEXP coerce_vec(SEXP x, SEXPTYPE type){
  return Rf_coerceVector(x, type);
}

// One-parameter template version
template <RType T>
inline SEXP new_vector_impl(R_xlen_t n) {
  static_assert(
    always_false<T>,
    "Unimplemented `new_vector_impl` specialisation"
  );
  return r_null;
}

template <>
inline SEXP new_vector_impl<r_bool_t>(R_xlen_t n) {
  return internal::new_vec(LGLSXP, n);
}
template <>
inline SEXP new_vector_impl<r_int_t>(R_xlen_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline SEXP new_vector_impl<r_double_t>(R_xlen_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline SEXP new_vector_impl<r_int64_t>(R_xlen_t n){
  SEXP out = SHIELD(internal::new_vec(REALSXP, n));
  attr::set_old_class(out, SHIELD(internal::make_utf8_strsxp("integer64")));
  YIELD(2);
  return out;
}
template <>
inline SEXP new_vector_impl<r_string_t>(R_xlen_t n){
  return internal::new_vec(STRSXP, n);
}
template <>
inline SEXP new_vector_impl<r_complex_t>(R_xlen_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline SEXP new_vector_impl<r_byte_t>(R_xlen_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline SEXP new_vector_impl<sexp_t>(R_xlen_t n){
  return internal::new_vec(VECSXP, n);
}

}

template<RType T>
struct r_vector_t {
  SEXP value;
  const T* const_ptr;  // Always created
  T* ptr;              // Only initialized if writable

  // Constructor that wraps new_vector<T>
  explicit r_vector_t(R_xlen_t size)
    : value(internal::new_vector_impl<std::remove_cvref_t<T>>(size))
    , const_ptr(vector_ptr<const T>(value))
    , ptr(nullptr)
  {
    if constexpr (is_r_ptr_writable_v<T>) {
      ptr = vector_ptr<T>(value);
    }
  }

  // Constructor from existing SEXP
  explicit r_vector_t(SEXP s)
    : value(s)
    , const_ptr(vector_ptr<const T>(s))
    , ptr(nullptr)
  {
    if constexpr (is_r_ptr_writable_v<T>) {
      ptr = vector_ptr<T>(s);
    }
  }

  // get uses const pointer
  T get(R_xlen_t index) const {
    return const_ptr[index];
  }

  // Set only available if writable
  void set(R_xlen_t index, T val) {
    if constexpr (is_r_ptr_writable_v<T>){
      ptr[index] = val;
    } else {
      internal::set_value<T>(value, index, val);
    }
  }

  // Direct pointer access
  T* data() requires is_r_ptr_writable_v<T>{
    return ptr;
  }

  const T* data() const {
    return const_ptr;
  }

  // Implicit conversion to SEXP
  operator SEXP() const {
    return value;
  }

  // Get size
  R_xlen_t size() const {
    return Rf_xlength(value);
  }

  void fill(R_xlen_t start, R_xlen_t n, const T val){
  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    auto *p_target = data();
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        p_target[start + i] = val;
      }
    } else {
      std::fill_n(p_target + start, n, val);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      set(start + i, val);
    }
  }
}

void replace(R_xlen_t start, R_xlen_t n, const T old_val, const T new_val){
  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    auto *p_target = data();
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_target[start + i] == old_val){
          p_target[start + i] = new_val;
        }
      }
    } else {
      std::replace(data() + start, data() + start + n, old_val, new_val);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      R_xlen_t idx = start + i;
      if (get(idx) == old_val){
        set(idx, new_val);
      }
    }
  }
}

};

template<typename T>
struct is_r_vector : std::false_type {};

template<typename T>
struct is_r_vector<r_vector_t<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_r_vector_v = is_r_vector<std::remove_cvref_t<T>>::value;

}

#endif
