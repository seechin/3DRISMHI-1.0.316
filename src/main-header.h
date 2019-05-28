
#define     _GROMACS4_
//#define     _FFTWMPPARALLEL_                // fftw multithreading
#define     _LOCALPARALLEL_                 // general local paralleling
#define     _LIBZ_                          // use libz to compress TENSOR4D
#define     _INTERACTIVE_                   // allow interactive mode, undef this will disable main-interactive.cpp

#if defined(_LOCALPARALLEL_) || !defined(_FFTWMPPARALLEL_)
  #define   _LOCALPARALLEL_FFTW_            // local paralleling with manually assign fftw to multithreads
#endif
#ifdef _LOCALPARALLEL_
  #define   _LOCALPARALLEL_PTHREAD_         // enable parallel by thread
#endif

//#undef      _TTYPROMPTCOLOR_
