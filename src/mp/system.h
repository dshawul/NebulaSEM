#ifndef __SYSTEM_H
#define __SYSTEM_H

#include <string>
#include <cstdarg>
#ifdef _MSC_VER
#    include <windows.h>
#    include <process.h>
#	 include <sys/timeb.h>
#else
#    include <unistd.h>
#    include <sys/stat.h>
#	 include <sys/time.h>
#endif


namespace System {
	/*get processor id*/
	inline int get_pid() {
#ifdef _MSC_VER
		return _getpid();
#else
		return getpid();
#endif
	}
	/*system dependent directory operations*/
	inline int cd(std::string path) {
#ifdef _MSC_VER
		return ::SetCurrentDirectory((LPCTSTR)path.c_str());
#else
		return !::chdir(path.c_str());
#endif
	}
	inline int mkdir(std::string path) {
#ifdef _MSC_VER
		return ::CreateDirectory((LPCTSTR)path.c_str(),NULL);
#else
		return !::mkdir(path.c_str(),S_IRWXU);
#endif
	}
	inline int rmdir(std::string path) {
#ifdef _MSC_VER
		return ::RemoveDirectory((LPCTSTR)path.c_str());
#else
		return !::rmdir(path.c_str());
#endif
	}
	/*time*/
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
