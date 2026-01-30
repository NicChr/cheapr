#ifndef CHEAPR_REGISTER_HPP
#define CHEAPR_REGISTER_HPP

// The R parser will search for the string "[[cpp20::register]]"
#ifdef __R_GENERATE_
  #define CHEAPR_REGISTER [[cpp20::register]]
#else
  #define CHEAPR_REGISTER 
#endif

#endif
