#ifndef __TYPES_H
#define __TYPES_H

#ifdef _MSC_VER
#    pragma warning (disable: 4996)
#endif

/** Unsigned integer type*/
typedef unsigned int  Int;

/** Forces inlining of fucntions*/
#ifdef _MSC_VER
#   define FORCEINLINE __forceinline
#else
#   define FORCEINLINE __inline
#endif

/*
cache line memory alignment (64 bytes)
*/
#include <cstdlib>
#define CACHE_LINE_SIZE  64
#define CACHE_ALIGN alignas(CACHE_LINE_SIZE)

template<typename T>
void aligned_free(T*& mem) {
#ifdef _WIN32
    if(mem) _aligned_free(mem);
#else
    if(mem) free(mem);
#endif
    mem = 0;
}

template<typename T, int ALIGNMENT = CACHE_LINE_SIZE, bool large_pages = false>
void aligned_reserve(T*& mem,const size_t& size) {
#ifdef __ANDROID__
    mem = (T*) memalign(ALIGNMENT,size * sizeof(T));
#elif defined(_WIN32)
    mem = (T*)_aligned_malloc(size * sizeof(T),ALIGNMENT);
#else
    int ret = posix_memalign((void**)&mem,ALIGNMENT,size * sizeof(T));
    (void) ret;
#if defined(MADV_HUGEPAGE)
    if(large_pages)
        madvise(mem,size * sizeof(T),MADV_HUGEPAGE);
#endif
#endif
}

#endif
