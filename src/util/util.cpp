#include "util.h"

using namespace std;

namespace Util {
    std::map<std::string,ParamList*> ParamList::list;
}

/** String hash function */
Int Util::hash_function(std::string s) {
    Int h = 0;
    const char* p = s.c_str();
    while(*p) { h = 31 * h + *p++; }
    return h;
}

/** Read next character */
int Util::nextc(std::istream& is) {
    char c;
    is >> c;
    while(c == '#') {
        while((c = is.get()) && c != '\n');
        is >> c;
    }
    if(is.eof())
        return 0;
    is.putback(c);
    return c;
}

/** Reads parameters */
void Util::read_params(std::istream& is, bool output, std::string block) {
    string str;
    char c;

#define READ() {                \
    c = Util::nextc(is);        \
    if(!c) goto END;            \
    else if(c == '}') {         \
        is >> c;                \
        break;                  \
    } else is >> str;           \
}

    while(true) {
        READ();
        is >> c;

        auto it = ParamList::list.find(str);
        if((it == ParamList::list.end()) || 
                (!block.empty() && compare(str,block))) {
            int braces = 1;
            while((c = Util::nextc(is))) {
                is >> c;
                if(c == '{') braces++;
                else if(c == '}') {
                    braces--;
                    if(!braces) break;
                }
            }
            continue;
        }

        if(output) cout << str << "\n{" << endl;
        ParamList* params = it->second;
        while(true) {
            READ();
            if(output) cout << str << " = ";
            params->read(is,str,output);
            if(output) cout << endl;
        }
        if(output) cout << "}\n";
        if(!block.empty())
            break;
    }
END:
    is.clear();
    is.seekg(0,ios::beg);
}

