#define __REAL__  double
#define MACHINE_REASONABLE_ERROR 1e-12

#define PI  3.1415926535897932384626433832795
#define EE  2.7182818284590452353602874713527
#define COULCOOEF 138.9354846

#ifndef MAX_SOL
    #define MAX_SOL             100     // Max atom site number
#endif
#ifndef MAX_CMD_PARAMS
    #define MAX_CMD_PARAMS      MAX_SOL // Max parameter number for a command
#endif
#ifndef MAX_THREADS
    #define MAX_THREADS         500     // Max number of forks or threads
#endif
#ifndef MAX_DIIS
    #define MAX_DIIS            100     // Max DIIS steps
#endif
#define MAX_INCLUDE_RECURSIVE   100     // maximum include recursive levels


#define MAX_PATH                1024    // Maximum length of filename/path strings
#define MAX_NAME                64      // Maximum length of molecule/atom names
#define MAX_WORD                1000    // Maximum number of words in line analysis

#define MAX_RENAME_COUNT        10000   // Maximum number of renames
#define MAX_GVV_FILES           3

#ifdef VERSION
    #define DISTRIBUTE_VERSION  VERSION
#else
    #define DISTRIBUTE_VERSION  "1.0.316"
#endif
#ifndef DISTRIBUTE_VERSION
    const char * szLicence = "";
#else
  const char * szLicence = "\
  # 3DRISM-HI is free software for non-commercial use. You can use, modify \n\
  # or redistribute for non-commercial purposes under the terms of the GNU \n\
  # Lesser General Public License version 3\n\
";
#endif
