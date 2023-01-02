#ifndef __SYSTEM_H
#define __SYSTEM_H

#include <string>
#include <cstdarg>
#ifdef _MSC_VER
#    include <windows.h>
#    include <io.h>
#    include <process.h>
#    include <sys/timeb.h>
#else
#    include <unistd.h>
#    include <sys/stat.h>
#    include <sys/time.h>
#endif

/**
  System dependent collection of functions
 */
namespace System {
    /** Gets processor id*/
    inline int get_pid() {
#ifdef _MSC_VER
        return _getpid();
#else
        return getpid();
#endif
    }
    /** Change directory*/
    inline int cd(std::string path) {
#ifdef _MSC_VER
        return ::SetCurrentDirectory((LPCTSTR)path.c_str());
#else
        return !::chdir(path.c_str());
#endif
    }
    /** Makes new directory */
    inline int mkdir(std::string path) {
#ifdef _MSC_VER
        return ::CreateDirectory((LPCTSTR)path.c_str(),NULL);
#else
        return !::mkdir(path.c_str(),S_IRWXU);
#endif
    }
    /** Removes directory */
    inline int rmdir(std::string path) {
#ifdef _MSC_VER
        return ::RemoveDirectory((LPCTSTR)path.c_str());
#else
        return !::rmdir(path.c_str());
#endif
    }
    /** Gets present working directory */
    inline int pwd(char* path, int len) {
#ifdef _MSC_VER
        return ::GetCurrentDirectory(len,(LPTSTR)path);
#else
        return !::getcwd(path, len);
#endif
    }
    /** Check if file exists */
    inline int exists(std::string path) {
#ifdef _MSC_VER
        return (_access(path.c_str(), 0) == 0);
#else
        return (access(path.c_str(), F_OK) == 0);
#endif
    }
    /** Gets time in milli-seconds */
    inline int get_time() {
#ifdef _MSC_VER
        timeb tb;
        ftime(&tb);
        return int(tb.time * 1000 + tb.millitm);
#else
        timeval tb;
        gettimeofday(&tb, NULL);
        return int(tb.tv_sec * 1000 + tb.tv_usec / 1000);
#endif
    }
}

#endif
