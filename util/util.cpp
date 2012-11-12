#include <string>
#include "util.h"

using namespace std;

namespace Util {
	bool Terminated = false;
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
static void read(istream& is,string str,bool check) {
	using namespace Util;
	if(IntParams::read(str,is));
	else if(ScalarParams::read(str,is));
	else if(StringParams::read(str,is));
	else if(OptionParams::read(str,is));
	else if(VectorParams::read(str,is));
	else if(TensorParams::read(str,is));
	else if(STensorParams::read(str,is));
	else if(VerticesParams::read(str,is));
	else if(check) {
		cout << "UNKNOWN";
	}
}
void Util::read_params(istream& is) {
	cout << "Reading paramters:" << endl;
	string str;
	while(Util::nextc(is)) {
		is >> str;
		cout << str << " = ";
		read(is,str,true);
		cout << endl;
	}
	cout << "Finished reading." << endl;
}
void Util::read_param(istream& is,string param) {
	string str;
	while(Util::nextc(is)) {
		is >> str;
		if(!compare(str,param)) {
			cout << str << " = ";
			read(is,str,false);
			cout << endl;
			is.seekg(0,fstream::beg);
			break;
		}
	}
}
/*IntVector output*/
std::ostream& operator << (std::ostream& os, const std::vector<Int>& p) {
	Int sz = p.size();
	if(sz >= 16) os << sz << std::endl << "{ " << std::endl;
	else os << sz << "{ ";
	for(Int i = 0;i < sz;i++) {
		os << p[i] << " ";
		if(sz >= 16 && ((i + 1) % 16) == 0)
			os << std::endl;
	}
	if(sz >= 16) os << std::endl << "}";
	else os << "}";
	return os;
}
