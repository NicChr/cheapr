#ifndef CHEAPR_R_VECTOR_H
#define CHEAPR_R_VECTOR_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_attrs.h>
#include <cheapr/internal/r_vector_utils.h>
#include <cheapr/internal/r_rtype_coerce.h>

namespace cheapr {

template<RType T>
struct r_vec {
  SEXP value = R_NilValue;
  const T* const_ptr = nullptr;  // Always created
  T* ptr = nullptr;              // Only initialized if writable
  using data_type = T;           // Type of data vec contains

  // Constructor that wraps new_vector<T>
  explicit r_vec(r_size_t size)
    : value(internal::new_vector_impl<std::remove_cvref_t<T>>(size))
    , const_ptr(internal::vector_ptr<const T>(value))
    , ptr(nullptr)
  {
    if constexpr (RPtrWritableType<T>) {
      ptr = internal::vector_ptr<T>(value);
    }
  }

  // Constructor from existing SEXP
  explicit r_vec(SEXP s = R_NilValue)
    : value(s)
    , const_ptr(s == R_NilValue ? nullptr : internal::vector_ptr<const T>(s))
    , ptr(nullptr)
  {
    if constexpr (RPtrWritableType<T>){
      if (const_ptr != nullptr){
        ptr = internal::vector_ptr<T>(s);
      }
    }
  }

 // Maybe add implicit coercion to `r_sexp` in future
  // operator r_sexp() const {
  //   return r_sexp(value);
  // }
  // Implicit conversion to SEXP
  operator SEXP() const {
    return value;
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
    return value == R_NilValue;
  }

  // Get size
  r_size_t size() const {
    return Rf_xlength(value);
  }
  r_size_t length() const {
    return size();
  }

  const r_str address() const {
    return internal::address(value);
  }

  r_vec<r_str> get_names() const {
    return r_vec<r_str>(Rf_getAttrib(value, symbol::names_sym));
  }
  
  void set_names(r_vec<r_str> names){
    if (names.is_null()){
      attr::set_attr(value, symbol::names_sym, r_null);
    } else {
      Rf_namesgets(value, names);
    }
  }

  r_vec<r_sexp> get_attrs() const {
    return r_vec<r_sexp>(attr::get_attrs(value));
  }

  void clear_attrs(){
    auto attrs = SHIELD(get_attrs());
    if (attrs.is_null()){
      YIELD(1);
      return;
    }
    
    auto names = SHIELD(attrs.get_names());

    int n = attrs.length();
    for (int i = 0; i < n; ++i){
      r_sym target_sym = internal::as_r<r_sym>(names.get(i));
      attr::set_attr(value, target_sym, r_null);
    }
    YIELD(2);
}

  void set_attrs(r_vec<r_sexp> attrs) {

    if (this->is_null()){
      Rf_error("Cannot add attributes to `NULL`");
    }

    if (attrs.is_null()){
      return;
    }

    int32_t NP = 0;

    auto names = SHIELD(attrs.get_names()); ++NP;

    if (names.is_null()){
      YIELD(NP);
      Rf_error("attributes must be a named list");
    }

    clear_attrs();

    r_sym attr_nm;

    int n = names.length();

    for (int i = 0; i < n; ++i){
      if (!(names.get(i) == blank_r_string)){
        attr_nm = internal::as_r<r_sym>(names.get(i));
        if (this->address() == internal::address(static_cast<SEXP>(attrs.get(i)))){
          SEXP dup_attr = SHIELD(Rf_duplicate(attrs.get(i))); ++NP;
          attr::set_attr(value, attr_nm, dup_attr);
        } else {
          attr::set_attr(value, attr_nm, attrs.get(i));
        }
      }
    }
    YIELD(NP);
  }

  r_vec<r_lgl> is_na() const {
    r_size_t n = length();
    auto out = SHIELD(r_vec<r_lgl>(n));
    
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
    
    YIELD(1);
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
      internal::set_value<T>(value, index, val2);
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
      auto resized_vec = SHIELD(r_vec<T>(n));
      r_size_t n_to_copy = std::min(n, vec_size);
      std::copy_n(this->begin(), n_to_copy, resized_vec.begin());
      YIELD(1);
      return resized_vec;
    }
  }

  // Define this after make_vec<> and arg struct
  template<typename... Args>
  void modify_attrs(Args... args);

};


// R Date vector
struct r_dates : public r_vec<r_int> {
  
  // Constructors
  r_dates() : r_vec<r_int>() {}
  
  explicit r_dates(SEXP x) : r_vec<r_int>(x) {
    if (!is_null() && !attr::inherits1(x, "Date")){ 
      Rf_error("`SEXP` must be a Date");
    }
  }
  
  explicit r_dates(r_size_t n) : r_vec<r_int>(n) {
    attr::set_old_class(this->value, internal::make_utf8_strsxp("Date"));
  }
};

// R POSIXct vector
struct r_posixcts : public r_vec<r_dbl> {
  
  // Constructors
  r_posixcts() : r_vec<r_dbl>() {}
  
  explicit r_posixcts(SEXP x) : r_vec<r_dbl>(x) {
    if (!(is_null() || (attr::inherits1(x, "POSIXct") && attr::inherits1(x, "POSIXt")))){ 
      Rf_error("`SEXP` must be a POSIXct");
    }
  }
  
  explicit r_posixcts(r_size_t n) : r_vec<r_dbl>(n) {
    auto cls = SHIELD(r_vec<r_str>(2));
    auto tz = SHIELD(internal::new_scalar_vector(blank_r_string));
    cls.set(0, "POSIXct"); cls.set(1, "POSIXt");
    // Set class
    attr::set_old_class(this->value, cls);
    // Set timezone
    attr::set_attr(this->value, internal::as_r<r_sym>("tzone"), tz);
    YIELD(2);
  }

  r_str tzone_get() const {
    auto tz = r_vec<r_str>(attr::get_attr(this->value, internal::as_r<r_sym>("tzone")));
    if (tz.is_null() || tz.length() == 0){
      return blank_r_string;
    } else {
      return tz.get(0);
    }
  }

    void tzone_set(r_str tzone) {
    auto tz = SHIELD(r_vec<r_str>(attr::get_attr(this->value, internal::as_r<r_sym>("tzone"))));
    if (tz.length() != 0){
      tz.set(0, tzone);
      YIELD(1);
    } else {
      auto new_tz = SHIELD(internal::new_scalar_vector(tzone));
      attr::set_attr(this->value, internal::as_r<r_sym>("tzone"), new_tz);
      YIELD(2);
    }
  }

};


// R data frame
// struct r_df : public r_vec<r_sexp> {
  
//   // Constructors
//   r_df() : r_vec<r_sexp>() {}
  
//   explicit r_df(SEXP x) : r_vec<r_sexp>(x) {
//     if (!is_null() && !attr::inherits1(x, "Date")){ 
//       Rf_error("`SEXP` must be a Date");
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

}

template<typename T>
concept RVectorType = internal::is_r_vector_v<T> || is<T, r_dates> || is<T, r_posixcts>;

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

namespace vec {

// Templates for creating new vectors (can also be done via r_vec)
template <RType T>
inline r_vec<T> new_vector(r_size_t n){
  return r_vec<T>(n);
}

template <RType T, typename U>
inline r_vec<T> new_vector(r_size_t n, const U default_value){
  auto out = SHIELD(r_vec<T>(n));
  out.fill(0, n, default_value); 
  YIELD(1);
  return out;
}

inline bool is_object(SEXP x){
  return Rf_isObject(x);
}

inline bool is_bare(SEXP x){
  return !is_object(x);
}


template<typename T>
inline auto as_vector(T x){
  if constexpr (RVectorType<T>){
    return x;
  } else if constexpr (RType<T>){
    return new_vector<T>(1, x);
  } else if constexpr (IntegerType<T>){
    if constexpr (internal::can_be_int<T>){
      return new_vector<r_int>(1, r_int(static_cast<int>(x)));
    } else {
      return new_vector<r_dbl>(1, r_dbl(static_cast<double>(x)));
    }
  } else {
    static_assert(
      always_false<T>,
      "Unimplemented `as_vector` specialisation"
    );
  }
}
template<>
inline auto as_vector<bool>(bool x){
  return as_vector(r_lgl(x));
}
template<>
inline auto as_vector<int>(int x){
  return as_vector(r_int(x));
}
template<>
inline auto as_vector<double>(double x){
  return as_vector(r_dbl(x));
}
template<>
inline auto as_vector<const char *>(const char *x){
  return as_vector(r_str(internal::make_utf8_charsxp(x))); 
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
    return static_cast<SEXP>(new_vector<r_sexp>(1, r_sexp(x)));
  }
  }
}

}

}

#endif
