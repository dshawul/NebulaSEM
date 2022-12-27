#ifndef __UTIL_H
#define __UTIL_H

#include "tensor.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cstdarg>

/** \name Container iterators */
//@{
#define forEach(field,i)                                \
    for(Int i = 0;i < (field).size();i++)

#define forEachR(field,i)                               \
    for(Int i = (field).size();i-- > 0;)

#define forEachS(field,i,strt)                          \
    for(Int i = strt;i < (field).size();i++)

#define forEachSR(field,i,strt)                         \
    for(Int i = (field).size();i-- > strt;)

#define forEachIt(field,it)                        \
    for(auto it = (field).begin(); it != (field).end(); ++it)
//@}

/** Copy collection */            
#define copyColl(src,trg)                               \
    std::copy(std::begin(src),std::end(src),std::back_inserter(trg));

/** Erase by value */   
#define eraseValue(field,val)                           \
    (field).erase(std::remove((field).begin(),(field).end(),val),(field).end());

/** std::pair I/O*/
template <class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1,T2>& p) {
    os << p.first << " " << p.second;
    return os;
}

template <class T1, class T2>
std::istream& operator >> (std::istream& is, std::pair<T1,T2>& p) {
    is >> p.first >> p.second;
    return is;
}

/** std::vector I/O*/
template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& p) {
    os << p.size() << std::endl;
    os << "{" << std::endl;
    forEach(p,i)
        os << p[i] << std::endl;
    os << "}" << std::endl;
    return os;
}

template <class T>
std::istream& operator >> (std::istream& is, std::vector<T>& p) {
    Int size;
    char symbol;
    is >> size >> symbol;
    p.resize(size);
    forEach(p,i)
        is >> p[i];
    is >> symbol;
    return is;
}

template<>
std::ostream& operator << (std::ostream& os, const std::vector<Int>& p);

/** std::map I/O */
template <class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::map<T1,T2>& p) {
    os << p.size() << std::endl;
    os << "{" << std::endl;
    forEachIt(p,it)
        os << it->first << " " << it->second << std::endl;
    os << "}" << std::endl;
    return os;
}

template <class T1, class T2>
std::istream& operator >> (std::istream& is, std::map<T1,T2>& p) {
    Int size;
    char symbol;
    is >> size >> symbol;
    p.resize(size);
    forEachIt(p,it)
        is >> it->first >> it->second;
    is >> symbol;
    return is;
}

/** Test if two Int vectors are equal. O(n^2) complexity but arrays
  are small since it is used for edges/faces/cells*/
FORCEINLINE bool equalSet(std::vector<Int>& v1,std::vector<Int>& v2) {
    forEach(v1,i) {
        bool found = false;
        forEach(v2,j) {
            if(v1[i] == v2[j]) {
                found = true;
                break;
            }
        }
        if(!found)
            return false;
    }
    return true;
}

/** Erase indices from a vector. It assumes the indices are already sorted */
template<typename T>
void erase_indices(std::vector<T>& data, const std::vector<Int>& indicesToDelete) {
    if(indicesToDelete.size() == 0)
        return;

    std::vector<T> temp;
    temp.reserve(data.size() - indicesToDelete.size());

    auto itBlockBegin = data.begin();
    forEachIt(indicesToDelete,it) {
        auto itBlockEnd = data.begin() + *it;
        if(itBlockBegin != itBlockEnd)
            std::copy(itBlockBegin, itBlockEnd, std::back_inserter(temp));
        itBlockBegin = itBlockEnd + 1;
    }
    auto itBlockEnd = data.end();
    if(itBlockBegin != itBlockEnd) 
        std::copy(itBlockBegin, itBlockEnd, std::back_inserter(temp));

    data = temp;
}

/** Compare pair using second value */
template<template <typename> class P = std::less >
struct compare_pair_second {
    template<class T1, class T2> bool operator()(
            const std::pair<T1, T2>& left, 
            const std::pair<T1, T2>& right
            ) {
        return P<T2>()(left.second, right.second);
    }
};

/** Class for some utililty functions */
namespace Util {
    Int hash_function(std::string s);
    int nextc(std::istream&);
    void cleanup();

    /** Compare two strings case-insensitive */
    inline int compare(std::string& s1,std::string s2) {
        std::string t1 = s1,t2 = s2;
        std::transform(t1.begin(),t1.end(),t1.begin(),toupper);
        std::transform(t2.begin(),t2.end(),t2.begin(),toupper);
        return (t1 != t2);
    }

    /** General string option list*/
    struct Option {
        Int* val;
        std::vector<std::string> list;
        Option(void* v,Int N, ...) {
            val = (Int*)v;
            std::string str;
            list.assign(N,"");
            va_list ap;
            va_start(ap, N);
            for(Int i = 0;i < N;i++) {
                str = va_arg(ap,char*);
                list[i] = str;
            }
            va_end(ap);
        }
        Int getID(std::string str) {
            forEach(list,i) {
                if(!Util::compare(list[i],str)) 
                    return i;
            }
            std::cout << "Unknown parameter : " << str << std::endl;
            return 0;
        }
        friend std::istream& operator >> (std::istream& is, Option& p) {
            std::string str;
            is >> str;
            *(p.val) = p.getID(str);
            return is;
        }
        friend std::ostream& operator << (std::ostream& os, const Option& p) {
            os << p.list[*(p.val)];
            return os;
        }
    };

    /** Special boolean option as a YES/NO */
    struct BoolOption : public Option {
        BoolOption(void* v) :
            Option(v,2,"NO","YES")
        {
        }
    };

    /** List of parameters in a group */
    template <typename T> 
    class Parameters{
        private:
            std::map<std::string,T*> list;
        public:
            void enroll(std::string str,T* addr) {
                list[str] = addr;
            }
            bool read(std::string str,std::istream& is,bool out) {
                auto it = list.find(str);
                if(it != list.end()) {
                    is >> *(it->second);
                    if(out) std::cout << *(it->second);
                    return true;
                }
                return false;
            }
    };
    extern void read_params(std::istream&, bool print, std::string block = "");

    /** Parameter list*/
    struct ParamList {
        std::string name;
        static std::map<std::string,ParamList*> list;

        ParamList(std::string n) : name(n) {
            list[name] = this;
        }
        ~ParamList() {
            list.erase(name);
        }

#define addParam(T,N)                               \
    Parameters<T> params_##N;                       \
    void enroll(std::string str,T* addr) {          \
        params_##N.enroll(str,addr);                \
    }
        addParam(Int,Int);
        addParam(Scalar,Scalar);
        addParam(Vector,Vector);
        addParam(STensor,STensor);
        addParam(Tensor,Tensor);
        addParam(std::string,string);
        addParam(Option,Option);
        addParam(std::vector<Int>,vec_int);
        addParam(std::vector<std::string>,vec_string);
        addParam(std::vector<Scalar>,vec_scalar);
        addParam(std::vector<Vector>,vec_vector);
#undef addParam

        void read(std::istream& is,std::string str,bool out) {
#define readp(N)  params_##N.read(str,is,out)
            if(readp(Int));
            else if(readp(string));
            else if(readp(Option));
            else if(readp(Scalar));
            else if(readp(Vector));
            else if(readp(Tensor));
            else if(readp(STensor));
            else if(readp(vec_int));
            else if(readp(vec_scalar));
            else if(readp(vec_vector));
            else if(readp(vec_string));
            else if(out) {
                std::cout << "UNKNOWN";
            }
#undef readp
        }
        void read(std::istream& is) {
            read_params(is,false,name);
        }
    };
    /*end*/
}

#endif
