#include <cstring>
#include "mesh.h"
#include "hexMesh.h"
#include "mshMesh.h"

using namespace std;

/*boundary*/
struct Bdry {
	string name;
	IntVector index;
};
/*generate mesh*/
int main(int argc,char* argv[]) {
	using namespace Mesh;
	using namespace Util;
	Vertices keys;
	vector<Bdry> Bdrys;
	MergeObject bMerge;
	string str;
	string default_name;
	char* file_name = argv[1];
	bool import = false;
	char c;

	/*command line arguments*/
	for(int i = 1;i < argc;i++) {
		if(!strcmp(argv[i],"-i")) {
			i++;
			import = true;
			file_name = argv[i];
		}
	}

	/*input stream*/
	ifstream input(file_name);

	/*import*/
	if(import) {
		input >> hex;
		mshMesh(input,gMesh);
		input >> dec;
		cout << hex;
		cout << gVertices << endl;
		cout << gFacets << endl;
		cout << gCells << endl;
		for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it)
			cout << it->first << " " << it->second << endl;
		cout << dec;
		return 0;
	}

	/*read key points*/
	if(Util::nextc(input))
		input >> keys;

	while((c = Util::nextc(input)) != 0) {
		char symbol;
		if(isdigit(c)) {
			/*read indices to keys*/
			IntVector index;
			input >> index;

			Vertices v(index.size(),Vector(0));
			for(Int i = 0;i < v.size();i++) 
				v[i] = keys[index[i]];

			IntVector n;
			Int type;
			vector<Scalar> s(12,Scalar(1));
			vector<Int> t(12);

			input >> str;
			if(!compare(str,"linear")) {
				input >> n;
				type = LINEAR;
				t.assign(12,type);
			} else { 
				if(!compare(str,"geometric")) type = GEOMETRIC;
				else if(!compare(str,"wall")) type = WALL;
				else if(!compare(str,"mixed")) type = MIXED;
				else return 0;
			
				//read divisions
				vector<Scalar> ts(s);
				vector<Int> tt(t);

				Int sz;
				input >> n;
				input >> sz >> symbol;

				if(type == MIXED) {
					for(Int i = 0;i < sz ;i++) {
						input >> symbol;
						switch(symbol) {
							case 'l':
							case 'L':
								type = LINEAR;
								break;
							case 'g':
							case 'G':
								type = GEOMETRIC;
								break;
							case 'w':
							case 'W':
								type = WALL;
								break;
						}
						tt[i] = type;
						input >> ts[i];
					}
				} else {
					for(Int i = 0;i < sz ;i++) 
						input >> ts[i];
					tt.assign(12,type);
				}

				input >> symbol;

				//assign to each side
				Int r = 12 / sz;
				for(Int i = 0;i < sz;i++) {
					for(Int j = 0;j < r;j++) {
						if(i * r + j < 12) {
							s[i * r + j] = ts[i];
							t[i * r + j] = tt[i];
						}
					}
				}
			}

			//curved edges
			static const int sides[12][2] = {
				{0,1}, {3,2}, {7,6}, {4,5},
				{0,3}, {1,2}, {5,6}, {4,7},
				{0,4}, {1,5}, {2,6}, {3,7}
			};
			vector<Edge> edges(12);
			for(Int i = 0;i < 12;i++) {
				edges[i].v[0] = v[sides[i][0]];
				edges[i].v[1] = v[sides[i][1]];
			}

			if((c = Util::nextc(input)) && (c == 'e')) {
				Int sz,side,key;
				input >> str;
				if(!compare(str,"edges")) {
					input >> sz >> symbol;
					for(Int i = 0;i < sz;i++) {
						input >> str >> side >> key;	
						Edge& e = edges[side];
						e.v[2] = keys[key];
						if(!compare(str,"arc")) {
							e.type = ARC;
						} else if(!compare(str,"cosine")) {
							e.type = COSINE;
						} else if(!compare(str,"quad")) {
							e.type = QUAD;
						} else {
							e.type = NONE;
						}
					}
					input >> symbol;
				} else {
					Bdry b;
					b.name = str;
					while((c = Util::nextc(input)) && isdigit(c)) {
						input >> b.index;
						Bdrys.push_back(b);
					}
				}
			}

			//generate mesh
			MeshObject mo;
			hexMesh(&n[0],&s[0],&t[0],&v[0],&edges[0],mo);
			merge(gMesh,bMerge,mo);
		} else {
			/*read boundaries*/
			Bdry b;
			input >> b.name;
			if(b.name == "default") {
				input >> default_name;
			} else {
				while((c = Util::nextc(input)) && isdigit(c)) {
					input >> b.index;
					Bdrys.push_back(b);
				}
			}
		}
	}
	/*merge boundary & internals*/
	merge(gMesh,bMerge);

	/*entities*/
	cout << hex;
	cout << gVertices << endl;
	cout << gFacets << endl;
	cout << gCells << endl;
    
	/*boundaries*/
	Int i,j;
	for(i = 0;i < Bdrys.size();i++) {
		IntVector list;
		IntVector& b = Bdrys[i].index;
		Vector N = (keys[b[1]] - keys[b[0]]) ^ (keys[b[2]] - keys[b[0]]);
		N /= mag(N);
		for(j = gMesh.nf;j < gFacets.size();j++) {
			Facet& f = gFacets[j];
			Vector N1 = ((gVertices[f[1]] - gVertices[f[0]]) ^ (gVertices[f[2]] - gVertices[f[0]]));
			N1 /= mag(N1);
			Vector H = (gVertices[f[0]] - keys[b[0]]);
			Scalar d = mag(N ^ N1);
			Scalar d2 = sqrt(mag(N & H));
			if(d <= 10e-4 && d2 <= 10e-4) {
				list.push_back(j);
			}
		}
		if(!list.empty()) {
			IntVector& gB = gBoundaries[Bdrys[i].name.c_str()];
			IntVector::iterator it = find(gB.begin(),gB.end(),list[0]);
			if(it == gB.end()) {
				for(j = 0;j < list.size();j++)
					gB.push_back(list[j]);
			}
		}
	}
	/*default specified*/
	if(!default_name.empty()) {
		IntVector& gB = gBoundaries[default_name.c_str()];
		for(i = gMesh.nf;i < gFacets.size();i++) {
			if(!faceInBoundary(i)) {
				gB.push_back(i);
			}
		}
	}
	/*write boundaries*/
	for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it)
		cout << it->first << " " << it->second << endl;
	cout << dec;

	return 0;
}
