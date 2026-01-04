#ifndef CHEAPR_R_VECTOR_H
#define CHEAPR_R_VECTOR_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>

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
  SEXP out = Rf_protect(internal::new_vec(REALSXP, n));
  SEXP cls = Rf_protect(internal::make_utf8_strsxp("integer64"));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  Rf_unprotect(2);
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

namespace internal {

inline r_int_t* integer_ptr(SEXP x){
  return reinterpret_cast<r_int_t*>(INTEGER(x));
}
inline const r_int_t* integer_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int_t*>(INTEGER_RO(x));
}
inline r_bool_t* logical_ptr(SEXP x){
  return reinterpret_cast<r_bool_t*>(integer_ptr(x));
}
inline const r_bool_t* logical_ptr_ro(SEXP x){
  return reinterpret_cast<const r_bool_t*>(integer_ptr_ro(x));
}
inline r_double_t* real_ptr(SEXP x){
  return reinterpret_cast<r_double_t*>(REAL(x));
}
inline const r_double_t* real_ptr_ro(SEXP x){
  return reinterpret_cast<const r_double_t*>(REAL_RO(x));
}
inline r_int64_t* integer64_ptr(SEXP x){
  return reinterpret_cast<r_int64_t*>(real_ptr(x));
}
inline const r_int64_t* integer64_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int64_t*>(real_ptr_ro(x));
}
inline r_complex_t* complex_ptr(SEXP x){
  return reinterpret_cast<r_complex_t*>(COMPLEX(x));
}
inline const r_complex_t* complex_ptr_ro(SEXP x){
  return reinterpret_cast<const r_complex_t*>(COMPLEX_RO(x));
}
inline r_byte_t* raw_ptr(SEXP x){
  return reinterpret_cast<r_byte_t*>(RAW(x));
}
inline const r_byte_t* raw_ptr_ro(SEXP x){
  return reinterpret_cast<const r_byte_t*>(RAW_RO(x));
}
inline const r_string_t* string_ptr_ro(SEXP x){
  return reinterpret_cast<const r_string_t*>(STRING_PTR_RO(x));
}

}

template<RType T>
inline T* vector_ptr(SEXP x) {
  static_assert(
    always_false<T>,
    "Unsupported type for vector_ptr"
  );
  return nullptr;
}


template<>
inline r_bool_t* vector_ptr<r_bool_t>(SEXP x) {
  return reinterpret_cast<r_bool_t*>(LOGICAL(x));
}

template<>
inline const r_bool_t* vector_ptr<const r_bool_t>(SEXP x) {
  return reinterpret_cast<const r_bool_t*>(LOGICAL_RO(x));
}
template<>
inline r_int_t* vector_ptr<r_int_t>(SEXP x) {
  return reinterpret_cast<r_int_t*>(INTEGER(x));
}
template<>
inline const r_int_t* vector_ptr<const r_int_t>(SEXP x) {
  return reinterpret_cast<const r_int_t*>(INTEGER_RO(x));
}

template<>
inline r_double_t* vector_ptr<r_double_t>(SEXP x) {
  return reinterpret_cast<r_double_t*>(REAL(x));
}
template<>
inline const r_double_t* vector_ptr<const r_double_t>(SEXP x) {
  return reinterpret_cast<const r_double_t*>(REAL_RO(x));
}

template<>
inline r_int64_t* vector_ptr<r_int64_t>(SEXP x) {
  return reinterpret_cast<r_int64_t*>(REAL(x));
}
template<>
inline const r_int64_t* vector_ptr<const r_int64_t>(SEXP x) {
  return reinterpret_cast<const r_int64_t*>(REAL_RO(x));
}


template<>
inline r_complex_t* vector_ptr<r_complex_t>(SEXP x) {
  return reinterpret_cast<r_complex_t*>(COMPLEX(x));
}
template<>
inline const r_complex_t* vector_ptr<const r_complex_t>(SEXP x) {
  return reinterpret_cast<const r_complex_t*>(COMPLEX_RO(x));
}

template<>
inline r_byte_t* vector_ptr<r_byte_t>(SEXP x) {
  return reinterpret_cast<r_byte_t*>(RAW(x));
}
template<>
inline const r_byte_t* vector_ptr<const r_byte_t>(SEXP x) {
  return reinterpret_cast<const r_byte_t*>(RAW_RO(x));
}

template<>
inline const r_string_t* vector_ptr<const r_string_t>(SEXP x) {
  return reinterpret_cast<const r_string_t*>(STRING_PTR_RO(x));
}

template<>
inline const sexp_t* vector_ptr<const sexp_t>(SEXP x) {
  return reinterpret_cast<const sexp_t*>(VECTOR_PTR_RO(x));
}

namespace internal {


// R vector getters + setters

template<RType T>
inline T get_value(SEXP x, const R_xlen_t i){
  return vector_ptr<T>(x)[i];
}

template<RType T>
inline T get_value(T *p_x, const R_xlen_t i){
  return p_x[i];
}

template<RType T>
inline const T get_value(const T *p_x, const R_xlen_t i){
  return p_x[i];
}


template<typename T>
requires (RType<T> || is<T, const char *>)
inline void set_value(T *p_x, const R_xlen_t i, T val){
  p_x[i] = val;
}

template<typename T>
requires (RType<T> || is<T, const char *>)
inline void set_value(SEXP x, const R_xlen_t i, T val){
  vector_ptr<T>(x)[i] = val;
}
template<>
inline void set_value<r_complex_t>(r_complex_t* p_x, const R_xlen_t i, r_complex_t val){
  p_x[i].re() = val.re();
  p_x[i].im() = val.im();
}
template<>
inline void set_value<r_complex_t>(SEXP x, const R_xlen_t i, r_complex_t val){
  auto *p_x = vector_ptr<r_complex_t>(x);
  p_x[i].re() = val.re(); 
  p_x[i].im() = val.im();
}
template<>
inline void set_value<r_string_t>(SEXP x, const R_xlen_t i, r_string_t val){
  SET_STRING_ELT(x, i, static_cast<SEXP>(val));
}
template<>
inline void set_value<const char *>(SEXP x, const R_xlen_t i, const char* val){
  set_value<r_string_t>(x, i, static_cast<r_string_t>(internal::make_utf8_charsxp(val)));
}

// Never use the pointer here to assign
template<>
inline void set_value<sexp_t>(SEXP x, const R_xlen_t i, sexp_t val){
  SET_VECTOR_ELT(x, i, val);
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
  R_xlen_t length() const {
    return Rf_xlength(value);
  }

  const r_string_t address() const {
    return internal::address(value);
  }

  r_vector_t<r_bool_t> is_na() const {
    R_xlen_t n = length();
    auto out = SHIELD(r_vector_t<r_bool_t>(n));
    for (R_xlen_t i = 0; i < n; ++i){
      out.set(i, is_r_na(get(i)));
    }
    YIELD(1);
    return out;
  }

    // get uses const pointer
  T get(R_xlen_t index) const {
    return const_ptr[index];
  }

  // Set only available if writable
  // We use flexible template to be able to coerce it to an RType
  template <typename U>
   void set(R_xlen_t index, U val) {
    T val2 = as<T>(val);
    if constexpr (is_r_ptr_writable_v<T>){
      ptr[index] = val2;
    } else {
      internal::set_value<T>(value, index, val2);
    }
  }
  // void set(R_xlen_t index, T val) {
  //   if constexpr (is_r_ptr_writable_v<T>){
  //     ptr[index] = val;
  //   } else {
  //     internal::set_value<T>(value, index, val);
  //   }
  // }

  template <typename U>
  void fill(R_xlen_t start, R_xlen_t n, const U val){
  auto val2 = as<T>(val);
  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    auto *p_target = data();
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        p_target[start + i] = val2;
      }
    } else {
      std::fill_n(p_target + start, n, val2);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      set(start + i, val2);
    }
  }
}

template <typename U1, typename U2> 
void replace(R_xlen_t start, R_xlen_t n, const U1 old_val, const U2 new_val){
  auto old_val2 = as<T>(old_val);
  auto new_val2 = as<T>(new_val);
  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    auto *p_target = data();
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_target[start + i] == old_val2){
          p_target[start + i] = new_val2;
        }
      }
    } else {
      std::replace(data() + start, data() + start + n, old_val2, new_val2);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      R_xlen_t idx = start + i;
      if (get(idx) == old_val2){
        set(idx, new_val2);
      }
    }
  }
}

//   void fill(R_xlen_t start, R_xlen_t n, const T val){
//   if constexpr (is_r_ptr_writable_v<T>){
//     int n_threads = internal::calc_threads(n);
//     auto *p_target = data();
//     if (n_threads > 1) {
//       OMP_PARALLEL_FOR_SIMD(n_threads)
//       for (R_xlen_t i = 0; i < n; ++i) {
//         p_target[start + i] = val;
//       }
//     } else {
//       std::fill_n(p_target + start, n, val);
//     }
//   } else {
//     for (R_xlen_t i = 0; i < n; ++i) {
//       set(start + i, val);
//     }
//   }
// }

// void replace(R_xlen_t start, R_xlen_t n, const T old_val, const T new_val){
//   if constexpr (is_r_ptr_writable_v<T>){
//     int n_threads = internal::calc_threads(n);
//     auto *p_target = data();
//     if (n_threads > 1) {
//       OMP_PARALLEL_FOR_SIMD(n_threads)
//       for (R_xlen_t i = 0; i < n; ++i) {
//         if (p_target[start + i] == old_val){
//           p_target[start + i] = new_val;
//         }
//       }
//     } else {
//       std::replace(data() + start, data() + start + n, old_val, new_val);
//     }
//   } else {
//     for (R_xlen_t i = 0; i < n; ++i) {
//       R_xlen_t idx = start + i;
//       if (get(idx) == old_val){
//         set(idx, new_val);
//       }
//     }
//   }
// }

};

template<typename T>
struct is_r_vector : std::false_type {};

template<typename T>
struct is_r_vector<r_vector_t<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_r_vector_v = is_r_vector<std::remove_cvref_t<T>>::value;

namespace vec {

// Templates for creating new vectors (can also be done via r_vector_t)
template <RType T>
inline r_vector_t<T> new_vector(R_xlen_t n){
  return r_vector_t<T>(n);
}

template <RType T>
inline r_vector_t<T> new_vector(R_xlen_t n, const T default_value){
  auto out = SHIELD(r_vector_t<T>(n));
  out.fill(0, n, default_value); 
  YIELD(1);
  return out;
}

inline bool is_object(SEXP x){
  return Rf_isObject(x);
}

inline bool is_atomic(SEXP x){
  return Rf_isVectorAtomic(x);
}

inline bool is_vec(SEXP x){
  return Rf_isVector(x);
}

inline bool is_bare(SEXP x){
  return !is_object(x);
}

inline bool is_logical(SEXP x){
  return TYPEOF(x) == LGLSXP;
}
inline bool is_integer(SEXP x){
  return TYPEOF(x) == INTSXP;
}
inline bool is_double(SEXP x){
  return TYPEOF(x) == REALSXP;
}
inline bool is_integer64(SEXP x){
  return is_double(x) && attr::inherits1(x, "integer64");
}
inline bool is_character(SEXP x){
  return TYPEOF(x) == STRSXP;
}
inline bool is_list(SEXP x){
  return TYPEOF(x) == VECSXP;
}
inline bool is_complex(SEXP x){
  return TYPEOF(x) == CPLXSXP;
}
inline bool is_raw(SEXP x){
  return TYPEOF(x) == RAWSXP;
}
inline bool is_date(SEXP x){
  return attr::inherits1(x, "Date");
}
inline bool is_datetime(SEXP x){
  return attr::inherits1(x, "POSIXt");
}
inline bool is_factor(SEXP x){
  return is_integer(x) && attr::inherits1(x, "factor");
}


template<typename T>
inline auto as_vector(T x){
  if constexpr (RType<T>){
    return new_vector<T>(1, x);
  } else if constexpr (is_r_vector_v<T>){
    return x;
    // This is a catch-all for all C integer types that aren't int or int64_t which are handled via specialisations
    // Also they can't have NA values so it's okay to static cast
  } else if constexpr (is_r_or_cpp_integral_v<T>){
    if constexpr (internal::can_be_int<T>){
      return new_vector<r_int_t>(1, r_int_t(static_cast<int>(x)));
    } else {
      return new_vector<r_double_t>(1, r_double_t(static_cast<double>(x)));
    }
  } else {
    static_assert(
      always_false<T>,
      "Unimplemented `as_vector` specialisation"
    );
    return r_null;
  }
}
template<>
inline auto as_vector<bool>(bool x){
  return as_vector(r_bool_t(x));
}
template<>
inline auto as_vector<int>(int x){
  return as_vector(r_int_t(x));
}
template<>
inline auto as_vector<double>(double x){
  return as_vector(r_double_t(x));
}
template<>
inline auto as_vector<const char *>(const char *x){
  return as_vector(r_string_t(internal::make_utf8_charsxp(x)));
}

template<>
inline auto as_vector<SEXP>(SEXP x){ 
  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP:
  case REALSXP:
  case STRSXP:
  case VECSXP:
  case CPLXSXP:
  case RAWSXP: {
    return x;
  }
  default: {
    return static_cast<SEXP>(new_vector<sexp_t>(1, sexp_t(x)));
  }
  }
}

// Named argument

struct arg {
  const char* name;
  SEXP value;

  explicit arg(const char* n) : name(n), value(r_null) {}
  explicit arg(const char* n, SEXP v) : name(n), value(v) {}

  template<typename T>
  arg operator=(T v) const {
    return arg(name, as<SEXP>(v));
  }
};

// Variadic list constructor
template<typename... Args>
inline r_vector_t<sexp_t> make_list(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return r_vector_t<sexp_t>(0);
  } else {

    r_vector_t<sexp_t> out(n);

    // Are any args named?
    constexpr bool any_named = (is<Args, arg> || ...);

    r_vector_t<r_string_t> nms(any_named ? n : 0);
    SHIELD(nms);

    int i = 0;
    (([&]() {
      if constexpr (is<Args, arg>) { 
        out.set(i, args.value);
        nms.set(i, as<r_string_t>(args.name));
      } else {
        out.set(i, as<SEXP>(args));
      }
      ++i;
    }()), ...);

    if (any_named){
      attr::set_old_names(out, nms);
    }
    YIELD(2);
    return out;
  }
}

}

}

#endif
