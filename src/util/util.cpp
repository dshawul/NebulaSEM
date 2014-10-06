#include <string>
#include "util.h"

using namespace std;

namespace Util {
	bool Terminated = false;

	std::map<std::string,ParamList*> ParamList::list;
}

Int Util::hash_function(std::string s) {
	Int h = 0;
	const char* p = s.c_str();
	while(*p) { h = 31 * h + *p++; }
	return h;
}
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
void Util::cleanup () {
	Terminated = true;
	printf("Exiting application\n");
}
void Util::read_params(istream& is,std::string block) {
	string str;
	char c;
	bool output = block.empty();

#define READ() {				\
	c = Util::nextc(is);		\
	if(!c) goto END;			\
	else if(c == '}') {			\
		is >> c;				\
		break;					\
	} else is >> str;			\
}

	while(true) {
		READ();
		is >> c;

		map<string,ParamList*>::iterator it = ParamList::list.find(str);
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

		if(output) cout << str << "\n{\n" << endl;
		ParamList* params = it->second;
		while(true) {
			READ();
			if(output) cout << "\t" << str << " = ";
			params->read(is,str,output);
			if(output) cout << endl;
		}
		if(output) cout << "}\n" << endl;
		if(!block.empty())
			break;
	}
END:
	is.clear();
	is.seekg(0,ios::beg);
}

/*IntVector cout*/
template<>
std::ostream& operator << (std::ostream& os, const std::vector<Int>& p) {
	Int sz = p.size();
	if(sz >= 16) os << sz << std::endl << "{ ";
	else os << sz << "{ ";
	for(Int i = 0;i < sz;i++) {
		if(sz >= 16 && (i % 16) == 0)
			os << std::endl;
		os << p[i] << " ";
	}
	if(sz >= 16) os << std::endl << "}";
	else os << "}";
	return os;
}

