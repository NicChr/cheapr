#ifndef CHEAPR_R_SEQ_H
#define CHEAPR_R_SEQ_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_coerce.h>

namespace cheapr {

    // r_vec<r_sexp> sequences(r_vec<r_int> size, r_vec<r_int> from, r_vec<r_int> by){

    //     int32_t NP = 0;
      
    //     int size_n = size.length();
    //     int from_n = from.length();
    //     int by_n = by.length();
    //     if (size_n > 0 && (from_n <= 0 || by_n <= 0)){
    //       stop("from and by must both have length > 0");
    //     }
      
    //     double out_size = cpp_sum(size);
    //     double min_size = cpp_min(size);
    //     if (is_na(out_size)){
    //       Rf_error("size must not contain NA values");
    //     }
    //     if (min_size < 0){
    //       Rf_error("size must be a vector of non-negative integers");
    //     }
    //     R_xlen_t interrupt_counter = 0;
    //     int start, increment, seq_size;
    //     SEXP out = r_null;
      
      
    //     if (as_list){
      
    //       out = SHIELD(new_list(size_n)); ++NP;
    //       SEXP curr_seq;
      
    //       PROTECT_INDEX curr_seq_idx;
    //       R_ProtectWithIndex(curr_seq = r_null, &curr_seq_idx); ++NP;
      
    //       if (size_n > 0){
    //         const int *p_size = integer_ptr_ro(size);
    //         const int *p_from = integer_ptr_ro(from);
    //         const int *p_by = integer_ptr_ro(by);
    //         for (int i = 0, bi = 0, fi = 0; i < size_n;
    //           bi = (++bi == by_n) ? 0 : bi,
    //           fi = (++fi == from_n) ? 0 : fi,
    //           ++i){
    //           seq_size = p_size[i];
    //           R_Reprotect(curr_seq = vec::new_vector<int>(seq_size), curr_seq_idx);
    //           int* RESTRICT p_curr_seq = integer_ptr(curr_seq);
    //           start = p_from[fi];
    //           increment = p_by[bi];
    //           if (is_na(start)){
    //             YIELD(NP);
    //             Rf_error("from contains NA values");
    //           }
    //           if (is_na(increment)){
    //             YIELD(NP);
    //             Rf_error("by contains NA values");
    //           }
    //           for (int j = 0; j < seq_size; ++j, ++interrupt_counter, start += increment){
    //             if (interrupt_counter == 100000000){
    //               R_CheckUserInterrupt();
    //               interrupt_counter = 0;
    //             }
    //             p_curr_seq[j] = start;
    //           }
    //           SET_VECTOR_ELT(out, i, curr_seq);
    //         }
    //       }
      
    //     } else {
      
    //       R_xlen_t index = 0;
    //       out = SHIELD(vec::new_vector<int>(out_size)); ++NP;
    //       int* RESTRICT p_out = integer_ptr(out);
      
    //       if (size_n > 0){
    //         const int *p_size = integer_ptr_ro(size);
    //         const int *p_from = integer_ptr_ro(from);
    //         const int *p_by = integer_ptr_ro(by);
    //         for (int i = 0, bi = 0, fi = 0; i < size_n;
    //           bi = (++bi == by_n) ? 0 : bi,
    //           fi = (++fi == from_n) ? 0 : fi,
    //           ++i){
    //           seq_size = p_size[i];
    //           start = p_from[fi];
    //           increment = p_by[bi];
    //           if (is_na(start)){
    //             YIELD(NP);
    //             Rf_error("from contains NA values");
    //           }
    //           if (is_na(increment)){
    //             YIELD(NP);
    //             Rf_error("by contains NA values");
    //           }
    //           for (int j = 0; j < seq_size; ++j, ++index, ++interrupt_counter, start += increment){
    //             if (interrupt_counter == 100000000){
    //               R_CheckUserInterrupt();
    //               interrupt_counter = 0;
    //             }
    //             p_out[index] = start;
    //           }
    //         }
    //       }
    //     }
      
    //     YIELD(NP);
    //     return out;
    //   }

}

#endif
