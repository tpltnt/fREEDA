
/* __BEGIN_DECLS should be used at the beginning of your declarations,
   so that C++ compilers don't mangle their names.  Use __END_DECLS at
   the end of C declarations. */
/* #undef __BEGIN_DECLS */
/* #undef __END_DECLS */
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS
# define __END_DECLS
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif


/* Where the binary files go. */
#define BINARYDIR "/usr/local/bin/"

/* Where the data files go. */
#define BITMAPDIR "/usr/local/share/freeda-2.0/ifREEDA/bitmaps/"

/* Define if debug output should be supported. */
#define DEBUG 1

/* Where the documentation files go. */
#define DOCDIR "/usr/local/share/freeda-2.0/ifREEDA/docs/"

/* Define to 1 if you have the <ieeefp.h> header file. */
/* #undef HAVE_IEEEFP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Where the language files go. */
#define LANGUAGEDIR "/usr/local/share/freeda-2.0/ifREEDA/lang/"

/* Define if debug code should be suppressed. */
/* #undef NDEBUG */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "nkriplani@ncsu.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "ifREEDA"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "ifREEDA 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "ifREEDA"

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */
