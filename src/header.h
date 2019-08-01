
#define     _GROMACS4_
//#define     _FFTWMPPARALLEL_              // fftw multithreading
#define     _LOCALPARALLEL_                 // general local paralleling
#define     _LIBZ_                          // use libz to compress TENSOR4D
#define     _INTERACTIVE_                   // allow interactive mode, undef this will disable main-interactive.cpp

#if defined(_LOCALPARALLEL_) || !defined(_FFTWMPPARALLEL_)
  #define   _LOCALPARALLEL_FFTW_            // local paralleling with one fftw per thread/fork
#endif
#ifdef _LOCALPARALLEL_
  #define   _LOCALPARALLEL_PTHREAD_         // enable parallel by thread, otherwise by fork
#endif

//#undef      _TTYPROMPTCOLOR_


#define _LOCALPARALLEL_FFTW_
#define _FFTWMPPARALLEL_
