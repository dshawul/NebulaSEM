#ifndef __TYPES_H
#define __TYPES_H

#ifdef _MSC_VER
#    pragma warning (disable: 4996)
#endif

/*define double precision*/
#define __DOUBLE

/*unsigned int*/
typedef unsigned int  Int;

/*force inline*/
#ifdef _MSC_VER
#	define FORCEINLINE __forceinline
#else
#	define FORCEINLINE __inline
#endif

#endif
