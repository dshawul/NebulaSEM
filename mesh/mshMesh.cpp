#include "mshMesh.h"
#include <sstream>

using namespace std;

void mshMesh(std::istream& is,Mesh::MeshObject& mo) {
	char symbol,c;
	int id,ND,zone,findex,lindex,
		type,bctype,ftype,etype,
		node_start = 0,facet_start = 0;

	/*read id*/
	while((c = Util::nextc(is)) != 0) {
		is >> symbol >> id;
		switch(id) {
		case 0x0: 
			do{ is >> c; } while(c != ')');
			break;
		case 0x2:
			is >> ND >> symbol;
			break;
		case 0x10:
			is >> symbol >> zone;
			is >> findex >> lindex;
			is >> type >> ND;
			is >> symbol >> symbol;
			if(zone != 0) {
				Vertex v;
				for(int i = findex;i <= lindex;i++) {
					is >> v;
					mo.v.push_back(v);
				}
				is >> symbol >> symbol;
			} else {
				node_start = findex;
			}
			break;
		case 0x12:
			is >> symbol >> zone;
			is >> findex >> lindex;
			is >> type;
			if((c = Util::nextc(is)) == ')');
			else is >> etype;
			is >> symbol >> symbol;

			while(symbol == '(') {
				do{ is >> c; } while(c != ')');
				is >> symbol;
			}
			if(zone == 0) {
				mo.c.resize(lindex);
			}
			break;
		case 0x13:
			is >> symbol >> zone;
			is >> findex >> lindex;
			is >> bctype;
			if((c = Util::nextc(is)) == ')') ftype = bctype;
			else is >> ftype;
			is >> symbol >> symbol;

			if(zone != 0) {
				std::stringstream name;
				name << "zone" << zone;
				IntVector& gB = mo.bdry[name.str().c_str()];

				Facet f;
				int n,c0,c1,k;
				for(int i = findex;i <= lindex;i++) {
					f.clear();
					is >> n;
					for(int j = 0;j < n;j++) {
						is >> k;
						f.push_back(k - node_start);
					}
					mo.f.push_back(f);
					gB.push_back(i - facet_start);

					is >> c0 >> c1;
					if(c0 == 0) {
						mo.fo.push_back(Constants::MAX_INT);
					} else {
						mo.fo.push_back(c0 - 1);
					}
					if(c1 == 0) {
						mo.fn.push_back(Constants::MAX_INT);
					} else {
                        mo.fn.push_back(c1 - 1);
					}
				}
				is >> symbol >> symbol;
			} else {
				facet_start = findex;
			}
			break;
		default:
			is >> symbol >> zone;
			while(symbol == '(') {
				do{ is >> c; } while(c != ')');
				is >> symbol;
			}
			is >> symbol;
			break;
		}
	}

	/*add cells*/
	Int co,cn;
	for(Int i = 0;i < mo.f.size();i++) {
		co = mo.fo[i];
		cn = mo.fn[i];
		if(co != Constants::MAX_INT)
			mo.c[co].push_back(i);
		if(cn != Constants::MAX_INT)
			mo.c[cn].push_back(i);
	}
}
