# Initialization of autoconf
AC_INIT([rismhi3d], [1.0.315])
AC_CONFIG_AUX_DIR([config_aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

# Check for C compiler
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)

# External library path select
AC_ARG_WITH([fftw],[AS_HELP_STRING([--with-fftw=fftw_path], [define the installation path of FFTW3])],[],[])
AC_ARG_WITH([gmx],[AS_HELP_STRING([--with-gmx=gmx_path], [define the installation path of GROMACS])],[],[])

# Features
#AC_ARG_WITH([interactive],[AS_HELP_STRING([--with-interactive], [allow interactive mode of rismhi3d])], [],[])
#AC_ARG_WITH([mt-fftw],[AS_HELP_STRING([--with-mt-fftw], [allow multithread FFTW. FFTW needs to be compiled with multhread enabled for this option])], [],[])

# Linux libraries
AC_CHECK_HEADERS([pthread.h], [
  LDFLAGS="$LDFLAGS -lpthread"
],[],[])
AC_CHECK_HEADERS([zlib.h], [
  LDFLAGS="$LDFLAGS -lz"
],[],[])

# FFTW library: must exist
if (test "x$with_fftw" == "x") then
  AC_CHECK_HEADERS([fftw3.h], [],[AC_MSG_ERROR(FFTW3 not found)],[])
  EXTRA_LDFLAGS="$EXTRA_LDFLAGS -lfftw3"
else
  EXTRA_CFLAGS="$EXTRA_CFLAGS -I$with_fftw/include"
  EXTRA_LDFLAGS="$EXTRA_LDFLAGS -L$with_fftw/lib -lfftw3"
  AC_CHECK_HEADERS([$with_fftw/include/fftw3.h], [],[AC_MSG_ERROR(FFTW3 not found)],[])
fi

# GROMACS library: optional
GROMACS_VER=0
if (test "x$with_gmx" == "x") then
  if (test "$GROMACS_VER" == 0) then
    AC_CHECK_HEADERS([gromacs/xtcio.h], [GROMACS_VER=4],[],[])
    if (test "$GROMACS_VER" != 0) then
      EXTRA_LDFLAGS="$EXTRA_LDFLAGS -lgmx -ldl"
    fi
  fi
  if (test "$GROMACS_VER" == 0) then
    AC_CHECK_HEADERS([gromacs/fileio/xtcio.h], [GROMACS_VER=16],[],[])
    if (test "$GROMACS_VER" != 0) then
      EXTRA_LDFLAGS="$EXTRA_LDFLAGS -lgromacs -ldl"
    fi
  fi
else
  if (test "$GROMACS_VER" = 0) then
    AC_CHECK_HEADERS([$with_gmx/include/gromacs/xtcio.h], [GROMACS_VER=4],[],[])
    if (test "$GROMACS_VER" != 0) then
      EXTRA_CFLAGS="$EXTRA_CFLAGS -I$with_gmx/include"
      EXTRA_LDFLAGS="$EXTRA_LDFLAGS -L$with_gmx/lib -lgmx -ldl"
    fi
  fi
  if (test "$GROMACS_VER" == 0) then
    AC_CHECK_HEADERS([$with_gmx/include/gromacs/fileio/xtcio.h], GROMACS_VER=16[],[],[])
    if (test "$GROMACS_VER" != 0) then
      EXTRA_CFLAGS="$EXTRA_CFLAGS -I$with_gmx/include"
      EXTRA_LDFLAGS="$EXTRA_LDFLAGS -L$with_gmx/lib -lgromacs -ldl"
    fi
  fi
fi

AC_SUBST(EXTRA_CFLAGS)
AC_SUBST(EXTRA_LDFLAGS)

if (test "$GROMACS_VER" == 16) then
  AC_DEFINE([GROMACS16], [16], [read XTC with GROMACS 2016+])
else
  if (test "$GROMACS_VER" == 4) then
    AC_DEFINE([GROMACS4], [4], [read XTC with GROMACS 4])
  fi
fi

# header file config.h
AC_CONFIG_HEADERS([config.h:config.hin])

# run config.status and generate Makefile
AC_CONFIG_FILES([Makefile])

# Generate the output
AC_OUTPUT
