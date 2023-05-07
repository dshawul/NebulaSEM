#ifndef __TYPES_H
#define __TYPES_H

#ifdef _MSC_VER
#    pragma warning (disable: 4996)
#endif

/** Defines double precision */
#define __DOUBLE

/** Unsigned integer type*/
typedef unsigned int  Int;

/** Forces inlining of fucntions*/
#ifdef _MSC_VER
#   define FORCEINLINE __forceinline
#else
#   define FORCEINLINE __inline
#endif

#endif
