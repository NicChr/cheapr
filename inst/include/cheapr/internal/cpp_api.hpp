#ifndef CHEAPR_CPP_API_HPP
#define CHEAPR_CPP_API_HPP

#if defined(CHEAPR_API_SELECTED)
#if CHEAPR_API_SELECTED != 1
#error "Only one cheapr API header may be included"
#endif
#else
#define CHEAPR_API_SELECTED 1 // C++ API
#endif

#endif
