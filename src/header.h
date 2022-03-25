#include "../config.h"

// Linux libraries
#ifndef HAVE_PTHREAD_H 
  #define   _LOCALPARALLEL_                 // general local paralleling
#else
  #define   _LOCALPARALLEL_PTHREAD_
#endif
#ifdef HAVE_ZLIB_H
  #define   _LIBZ_                          // use libz to compress TENSOR4D
#endif

// Third party libraries
// Gromacs for XTC
#ifdef GROMACS4
  #define   _GROMACS4_
#endif
#ifdef GROMACS16
  #define   _GROMACS2016_
#endif

// Features
//#define     _INTERACTIVE_                 // allow interactive mode, undef this will disable main-interactive.cpp
#define     _TTYPROMPTCOLOR                 // allow color text in terminal_
//#define     _FUNCTION_EXPORT_
//#define     _FFTWMPPARALLEL_              // fftw multithreading
