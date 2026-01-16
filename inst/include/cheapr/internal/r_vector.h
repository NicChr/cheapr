#ifndef CHEAPR_R_VECTOR_H
#define CHEAPR_R_VECTOR_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_vector_utils.h>
#include <cheapr/internal/r_rtype_coerce.h>

namespace cheapr {

template<RType T>
struct r_vec {
  r_sexp sexp = r_null;
  const T* const_ptr = nullptr;  // Always created
  T* ptr = nullptr;              // Only initialized if writable
  using data_type = T;           // Type of data vec contains

  // Constructor that wraps new_vector_impl<T>
  explicit r_vec(r_size_t size)
    : sexp(internal::new_vector_impl<std::remove_cvref_t<T>>(size))
    , const_ptr(internal::vector_ptr<const T>(sexp))
    , ptr(nullptr)
  {
    // Rf_protect(sexp);
    if constexpr (RPtrWritableType<T>) {
      ptr = internal::vector_ptr<T>(sexp);
    }
  }

  template<typename U>
  explicit r_vec(r_size_t size, U default_value)
    : r_vec(size)
  {
    fill(0, size, default_value);
  }

  // Constructor from existing SEXP
  explicit r_vec(SEXP s = R_NilValue)
    : sexp(s)
    , const_ptr(s == R_NilValue ? nullptr : internal::vector_ptr<const T>(s))
    , ptr(nullptr)
  {
    // Rf_protect(sexp);
    if constexpr (RPtrWritableType<T>){
      if (const_ptr != nullptr){
        ptr = internal::vector_ptr<T>(s);
      }
    }
  }

  // // Destructor - automatic unprotection
  // ~r_vec() {
  //   Rf_unprotect(1);
  // }
  // // Copy constructor - shield the copy
  // r_vec(const r_vec& other)
  //   : sexp(other.sexp)
  //   , const_ptr(other.const_ptr)
  //   , ptr(other.ptr)
  // {
  //   Rf_protect(sexp);
  // }
  // // Copy assignment
  // r_vec& operator=(const r_vec& other) {
  //   if (this != &other) {
  //     Rf_unprotect(1);
  //     sexp = other.sexp;
  //     const_ptr = other.const_ptr;
  //     ptr = other.ptr;
  //     Rf_protect(sexp);
  //   }
  //   return *this;
  // }

  // // Move constructor - transfer ownership
  // r_vec(r_vec&& other) noexcept
  //   : sexp(other.sexp)
  //   , const_ptr(other.const_ptr)
  //   , ptr(other.ptr)
  // {
  //   other.sexp = r_null;
  //   other.const_ptr = nullptr;
  //   other.ptr = nullptr;
  // }

  // // Move assignment - transfer ownership
  // r_vec& operator=(r_vec&& other) noexcept {
  //   if (this != &other) {
  //     Rf_unprotect(1);
  //     sexp = other.sexp;
  //     const_ptr = other.const_ptr;
  //     ptr = other.ptr;
  //     other.sexp = r_null;
  //     other.const_ptr = nullptr;
  //     other.ptr = nullptr;
  //   }
  //   return *this;
  // }

  // Implicit conversion to SEXP
  constexpr operator SEXP() const {
    return sexp.value;
  }

  // Direct pointer access
  T* data() requires RPtrWritableType<T>{
    return ptr;
  }

  const T* data() const {
    return const_ptr;
  }

  // Iterator support - begin + end
  T* begin() requires RPtrWritableType<T> {
      return ptr;
  }

  T* end() requires RPtrWritableType<T> {
      return ptr + size();
  }

  const T* begin() const {
      return const_ptr;
  }

  const T* end() const {
      return const_ptr + size();
  }

  bool is_null() const {
    return sexp.is_null();
  }

  // Get size
  r_size_t size() const {
    return Rf_xlength(sexp);
  }
  r_size_t length() const {
    return size();
  }

  r_str address() const {
    return sexp.address();
  }

  bool is_bare(){
    return !Rf_isObject(sexp);
  }

  r_vec<r_str> names() const {
    return r_vec<r_str>(Rf_getAttrib(sexp, symbol::names_sym));
  }

  void set_names(r_vec<r_str> names){
    if (names.is_null()){
      Rf_setAttrib(sexp, symbol::names_sym, r_null);
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
          out.set(i, is_r_na(get(i)));
        }
      } else {
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
          out.set(i, is_r_na(get(i)));
        }
      }
    } else {
      for (r_size_t i = 0; i < n; ++i){
        out.set(i, is_r_na(get(i)));
      }
    }

    return out;
  }

    // get uses const pointer
  T get(r_size_t index) const {
    return const_ptr[index];
  }

  // Set only available if writable
  // We use flexible template to be able to coerce it to an RType
  template <typename U>
  void set(r_size_t index, U val) {
      T val2 = internal::as_r<T>(val);
      if constexpr (RPtrWritableType<T>){
      ptr[index] = val2;
      } else {
      internal::set_value<T>(sexp, index, val2);
      }
  }

  template <typename U>
  void fill(r_size_t start, r_size_t n, const U val){
    auto val2 = internal::as_r<T>(val);
    if constexpr (RPtrWritableType<T>){
      int n_threads = internal::calc_threads(n);
      auto *p_target = data();
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
    bool implicit_na_coercion = !is_r_na(old_val) && is_r_na(old_val2);
    if (!implicit_na_coercion){
      if constexpr (RPtrWritableType<T>){
        int n_threads = internal::calc_threads(n);
        auto *p_target = data();
        if (n_threads > 1) {
          OMP_PARALLEL_FOR_SIMD(n_threads)
          for (r_size_t i = 0; i < n; ++i) {
            if (p_target[start + i] == old_val2){
              p_target[start + i] = new_val2;
            }
          }
        } else {
          std::replace(data() + start, data() + start + n, old_val2, new_val2);
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
      std::copy_n(this->begin(), n_to_copy, resized_vec.begin());
      return resized_vec;
    }
  }

};


// R data frame
// struct r_df : public r_vec<r_sexp> {

//   // Constructors
//   r_df() : r_vec<r_sexp>() {}

//   explicit r_df(SEXP x) : r_vec<r_sexp>(x) {
//     if (!is_null() && !attr::inherits1(x, "Date")){
//       Rf_error("`SEXP` must be a data frame");
//     }
//   }

//   explicit r_df(r_size_t n) : r_vec<r_sexp>(n) {
//     attr::set_old_class(this->value, internal::make_utf8_strsxp("Date"));
//   }
// };

namespace internal {

template<typename T>
struct is_r_vector : std::false_type {};

template<typename T>
struct is_r_vector<r_vec<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_r_vector_v = is_r_vector<std::remove_cvref_t<T>>::value;

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

template <RType T>
inline void r_copy_n(r_vec<T> &target, r_vec<T> &source, r_size_t target_offset, r_size_t n){
  auto *p_source = source.data();

  if constexpr (RPtrWritableType<T>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (r_size_t i = 0; i < n; ++i) {
        target.set(target_offset + i, p_source[i]);
      }
    } else {
      std::copy_n(p_source, n, target.data() + target_offset);
    }
  } else {
    for (r_size_t i = 0; i < n; ++i) {
      target.set(target_offset + i, p_source[i]);
    }
  }
}

// template<typename T>
// inline auto as_vector(T x){
//   if constexpr (RVectorType<T>){
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
//     auto rt_val = internal::as_r_type(x);
//     return r_vec<decltype(rt_val)>(1, rt_val);
//   }
// }

}

#endif
