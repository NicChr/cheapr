#ifndef CHEAPR_R_VECTOR_H
#define CHEAPR_R_VECTOR_H

#include <cheapr/internal/r_methods.h>
#include <cheapr/internal/r_vec_utils.h>
#include <cheapr/internal/r_rtype_coerce.h>

namespace cheapr {

template<RVal T>
struct r_vec {

  private:

  // Initialise read-only ptr to: 
  // SEXP - If T is `r_sexp` or `r_str`
  // T - Otherwise
  using ptr_t = std::conditional_t<RPtrWritableType<T>, unwrapped_t<T>*, const SEXP*>;  
  ptr_t m_ptr = nullptr;

  public: 

  r_sexp sexp = r_null;

  using data_type = T;

  // Constructor that wraps new_vec_impl<T>
  explicit r_vec(r_size_t n)
    : sexp(internal::new_vec_impl<std::remove_cvref_t<T>>(n))
  {
    if constexpr (RPtrWritableType<T>) {
      m_ptr = internal::vector_ptr<T>(sexp.value);
    } else if constexpr (any<T, r_sexp, r_sym>) {
      m_ptr = (const SEXP*) VECTOR_PTR_RO(sexp.value);
      } else if constexpr (is<T, r_str>){
      m_ptr = (const SEXP*) STRING_PTR_RO(sexp.value);
      }
  }

  template<typename U>
  explicit r_vec(r_size_t n, U default_value)
    : r_vec(n)
  {
    fill(0, n, default_value);
  }
  
  r_vec(): r_vec(r_size_t(0)){}

  // Constructors from existing r_sexp/SEXP
  explicit r_vec(r_sexp s) : sexp(std::move(s)) {
    if (sexp.value != R_NilValue) {
      // vector_ptr helper must be updated to return SEXP* for r_sexp/r_str/r_sym
      // We cast strictly to the stored type
      if constexpr (RPtrWritableType<T>) {
        m_ptr = internal::vector_ptr<T>(sexp);
      } else if constexpr (any<T, r_sexp, r_sym>) {
        m_ptr = (const SEXP*) VECTOR_PTR_RO(sexp); 
      } else if constexpr (is<T, r_str>){
        m_ptr = (const SEXP*) STRING_PTR_RO(sexp);
      }
    }
}

explicit r_vec(SEXP s) : r_vec(r_sexp(s)) {}

  // Implicit conversion to SEXP
  operator SEXP() const {
    return sexp;
  }

  // Explicit conversion to r_sexp
  constexpr explicit operator r_sexp() const {
    return sexp;
  }

  // Direct pointer access
  ptr_t data() const {
    return m_ptr;
  }

  // Iterator support - begin + end
  ptr_t begin() {
      return data();
  }

  ptr_t end() {
      return data() + size();
  }

  bool is_null() const {
    return sexp.is_null();
  }

  r_size_t length() const noexcept {
    return sexp.length();
  }

  // Get size
  r_size_t size() const noexcept {
    return length();
  }

  r_str address() const {
    return sexp.address();
  }

  bool is_bare(){
    return !Rf_isObject(sexp);
  }

  // get uses const pointer
  T get(r_size_t index) const {
    return T(m_ptr[index]);
  }

  // Set only available if writable
  // We use flexible template to be able to coerce it to an RVal
  template <typename U>
  void set(r_size_t index, U val) {
      T val2 = internal::as_r<T>(val);
      if constexpr (any<T, r_sexp, r_sym>){
        SET_VECTOR_ELT(sexp, index, val2);
      } else if constexpr (is<T, r_str>){
        SET_STRING_ELT(sexp, index, val2);
      } else {
        m_ptr[index] = unwrap(val2);
      }
  }

  // T operator[](r_size_t i) const {
  //     return get(i);
  // }

  // r_vec<T> operator[](const r_vec<r_int>& indices) const {
  //   r_size_t n_out = indices.length();
  //   r_vec<T> out(n_out);

  // }

  r_vec<r_str> names() const {
    return r_vec<r_str>(Rf_getAttrib(sexp, symbol::names_sym));
  }

  void set_names(r_vec<r_str> names){
    if (names.is_null()){
      Rf_setAttrib(sexp, symbol::names_sym, r_null);
    } else if (names.length() != sexp.length()){
      cpp11::stop("`length(names)` must equal `length(x)`");
    } else {
      Rf_namesgets(sexp, names);
    }
  }

  r_vec<r_lgl> is_na() const {
    r_size_t n = length();
    auto out = r_vec<r_lgl>(n);

    if constexpr (RPtrWritableType<T>){
      int n_threads = internal::calc_threads(n);
      if (n_threads > 1){
        OMP_PARALLEL_FOR_SIMD(n_threads)
        for (r_size_t i = 0; i < n; ++i){
          out.set(i, cheapr::is_na(get(i)));
        }
      } else {
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
          out.set(i, cheapr::is_na(get(i)));
        }
      }
    } else {
      for (r_size_t i = 0; i < n; ++i){
        out.set(i, cheapr::is_na(get(i)));
      }
    }

    return out;
  }

  template <typename U>
  void fill(r_size_t start, r_size_t n, const U val){
    auto val2 = unwrap(internal::as_r<T>(val));
    if constexpr (RPtrWritableType<T>){
      int n_threads = internal::calc_threads(n);
      auto* RESTRICT p_target = data();
      if (n_threads > 1) {
        OMP_PARALLEL_FOR_SIMD(n_threads)
        for (r_size_t i = 0; i < n; ++i) {
          p_target[start + i] = val2;
        }
      } else {
        std::fill_n(p_target + start, n, val2);
      }
    } else {
      for (r_size_t i = 0; i < n; ++i) {
        set(start + i, val2);
      }
    }
  }

  template <typename U1, typename U2>
  void replace(r_size_t start, r_size_t n, const U1 old_val, const U2 new_val){
    auto old_val2 = internal::as_r<T>(old_val);
    auto new_val2 = internal::as_r<T>(new_val);
    bool implicit_na_coercion = !cheapr::is_na(old_val) && cheapr::is_na(old_val2);
    if (!implicit_na_coercion){
      if constexpr (RPtrWritableType<T>){
        int n_threads = internal::calc_threads(n);
        auto *p_target = data();
        if (n_threads > 1) {
          OMP_PARALLEL_FOR_SIMD(n_threads)
          for (r_size_t i = 0; i < n; ++i) {
            r_lgl eq = p_target[start + i] == old_val2;
            if (eq.is_true()){
              p_target[start + i] = unwrap(new_val2);
            }
          }
        } else {
          OMP_SIMD
          for (r_size_t i = 0; i < n; ++i) {
            r_lgl eq = p_target[start + i] == old_val2;
            if (eq.is_true()){
              p_target[start + i] = unwrap(new_val2);
            }
          }
        }
      } else {
        for (r_size_t i = 0; i < n; ++i) {
          r_size_t idx = start + i;
          if (get(idx) == old_val2){
            set(idx, new_val2);
          }
        }
      }
    }
  }

  r_vec<T> resize(r_size_t n){
    r_size_t vec_size = length();
    if (n == vec_size){
      return *this;
    } else {
      auto resized_vec = r_vec<T>(n);
      r_size_t n_to_copy = std::min(n, vec_size);

      if constexpr (RPtrWritableType<T>){
        std::copy_n(this->begin(), n_to_copy, resized_vec.begin());
      } else {
        for (r_size_t i = 0; i < n_to_copy; ++i){
          resized_vec.set(i, this->get(i)); 
        }
      }
      return resized_vec;
    }
  }

  r_vec<T> rep_len(r_size_t n){

    r_size_t size = length();

    if (size == n){
      return *this;
    }

    auto out = r_vec<T>(n);

    if (size == 1){
      auto fill_val = this->get(0);
      out.fill(0, n, fill_val);
    } else if (n > 0 && size > 0){
      // Copy first block
      r_copy_n(out, *this, 0, size);

      // copy result to itself, doubling each iteration
      r_size_t copied = size;
      while (copied < n) {
        r_size_t to_copy = std::min(copied, n - copied);
        r_copy_n(out, out, copied, to_copy);
        copied += to_copy;
      }
      // If length > 0 but length(x) == 0 then fill with NA
    } else if (size == 0 && n > 0){
      auto na_val = na_value<T>();
      out.fill(0, n, na_val);
    }
    return out;
  }

};

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
  default:              cpp11::stop("`x` must be a vector");
  }
}

// Same as above but no run-time error, user must deal with non-vector input
template <class F>
decltype(auto) visit_maybe_vector(SEXP x, F&& f) {
  switch (CHEAPR_TYPEOF(x)) {
  case LGLSXP:          return f(r_vec<r_lgl>(x));
  case INTSXP:          return f(r_vec<r_int>(x));
  case CHEAPR_INT64SXP: return f(r_vec<r_int64>(x));
  case REALSXP:         return f(r_vec<r_dbl>(x));
  case STRSXP:          return f(r_vec<r_str>(x));
  case VECSXP:          return f(r_vec<r_sexp>(x));
  case CPLXSXP:         return f(r_vec<r_cplx>(x));
  case RAWSXP:          return f(r_vec<r_raw>(x));
  default:              return f(nullptr);
  }
}

}


template <RVal T>
inline void r_copy_n(r_vec<T> &target, r_vec<T> &source, r_size_t target_offset, r_size_t n){

  if constexpr (RPtrWritableType<T>){
    auto *p_source = source.data();
    auto *p_target = target.data();

    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (r_size_t i = 0; i < n; ++i) {
        p_target[target_offset + i] = p_source[i];
      }
    } else {
      std::copy_n(p_source, n, p_target + target_offset);
    }
  } else {
    for (r_size_t i = 0; i < n; ++i) {
      target.set(target_offset + i, source.get(i));
    }
  }
}

// template<typename T>
// inline auto as_vector(T x){
//   if constexpr (RVector<T>){
//     return x;
//   } else if constexpr (any<T, SEXP, r_sexp>){
//     switch (TYPEOF(x)){
//       case LGLSXP:
//       case INTSXP:
//       case REALSXP:
//       case STRSXP:
//       case VECSXP:
//       case CPLXSXP:
//       case RAWSXP: {
//         return r_sexp(x);
//       }
//       default: {
//         // New list of length 1 containing x
//         return r_sexp(static_cast<SEXP>(r_vec<r_sexp>(1, r_sexp(static_cast<SEXP>(x)))));
//       }
//       }
//   } else {
//     auto rt_val = internal::as_r_scalar(x);
//     return r_vec<decltype(rt_val)>(1, rt_val);
//   }
// }

}

#endif
