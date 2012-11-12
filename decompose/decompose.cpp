#include <sstream>
#include "decompose.h"
#include "system.h"

using namespace std;
using namespace Mesh;

/*duplicate file*/
template <class T>
void duplicate(istream& is,ostream& of) {
	MeshField<T,CELL> f;
	string str;
	int size;

	/*index file*/
	IntVector cLoc;
	ifstream index("index");
	index >> cLoc;

	/*read size*/
	is >> str >> size;
	of << str << " "<< size << endl;

	/*internal field*/
	char c;
	T value;
	if((c = Util::nextc(is)) && isalpha(c)) {
		value = T(0);
		is >> str;
		if(str == "uniform")
			is >> value;
		f = value;
	} else {
		char symbol;
		is >> size >> symbol;
		for(int i = 0;i < size;i++)
			is >> f[i];
		is >> symbol;
	}

	/*write it out*/
	of << cLoc.size() << endl;
	of << "{" << endl;
	for(Int j = 0;j < cLoc.size();j++)
		of << f[cLoc[j]] << endl;
	of << "}" << endl;

	/*boundaries*/
	string bname,cname;
	while((c = Util::nextc(is)) && isalpha(c)) {
		BCondition<T> bc(" ");
		is >> bc;
		of << bc << endl;
	}
	/*interMesh boundaries*/
	while((c = Util::nextc(index)) && isalpha(c)) {
		BCondition<T> bc(" ");
		index >> bc;
		of << bc << endl;
	}
}
/*decompose fields*/
void Decompose::decomposeFields(vector<string>& fields,std::string mName) {

	/*write meshes to file */
	int size;
	std::string str;

	for(Int ID = 0;;ID++) {
		/*cd*/
		stringstream path;
		path << mName << ID;
		if(!System::cd(path.str()))
			break;
		
		/*for each field*/
		for(Int i = 0;i < fields.size();i++) {
			/*read at time 0*/

			string str = "../" + fields[i] + "0";
			ifstream is(str.c_str());
			if(!is.fail()) {
				str = fields[i] + "0";
				ofstream of(str.c_str());

				/*seekg to beg*/
				is >> str >> size;
				is.seekg(0,fstream::beg);

				/*fields*/
				switch(size) {
				case 1 :  duplicate<Scalar>(is,of); break;
				case 3 :  duplicate<Vector>(is,of); break;
				case 6 :  duplicate<STensor>(is,of); break;
				case 9 :  duplicate<Tensor>(is,of); break;
				}
				/*end*/
			}
		}
		/*go back*/
		System::cd("..");
	}
}
/*decompose in x,y,z direction*/
int Decompose::decomposeXYZ(Mesh::MeshObject& mo,Int* n) {
	
	using Constants::MAX_INT;
	Int i,j,k,ID,count,total = n[0] * n[1] * n[2];
	Vector maxV(Scalar(-10e30)),minV(Scalar(10e30)),delta;

	/*facet owner and neighbors*/
	gFO.assign(mo.f.size(),MAX_INT);
	gFN.assign(mo.f.size(),MAX_INT);
	for(i = 0;i < mo.c.size();i++) {
		for(j = 0;j < mo.c[i].size();j++) {
			ID = mo.c[i][j];
			if(gFO[ID] == MAX_INT) 
				gFO[ID] = i;
			else 
				gFN[ID] = i;
		}
	}

	/*decomposed mesh*/
	MeshObject* meshes = new MeshObject[total];
	IntVector* vLoc = new IntVector[total];
	IntVector* fLoc = new IntVector[total];
	IntVector* cLoc = new IntVector[total];
	for(i = 0;i < total;i++) {
		vLoc[i].assign(mo.v.size(),0);
		fLoc[i].assign(mo.f.size(),0);
	}

	/*max and min points*/
	for(i = 0;i < mo.v.size();i++) {
		for(j = 0;j < 3;j++) {
			if(mo.v[i][j] > maxV[j]) maxV[j] = mo.v[i][j];
			if(mo.v[i][j] < minV[j]) minV[j] = mo.v[i][j];
		}
	}
    delta = maxV - minV;
	for(j = 0;j < 3;j++) delta[j] /= Scalar(n[j]);

	/*decompose cells*/
	MeshObject *pmesh;
	IntVector *pvLoc,*pfLoc,blockIndex;
	Vector C;
	
	blockIndex.assign(mo.c.size(),0);

	for(i = 0;i < mo.c.size();i++) {
		Cell& c = mo.c[i];

		/* approximate cell centre */
		C = Vector(0);
		count = 0;
		for(j = 0;j < c.size();j++) {
			Facet& f = mo.f[c[j]];
			for(k = 0;k < f.size();k++) {
				C += mo.v[f[k]];
				count++;
			}
		}
		C /= Scalar(count);

		/* add cell */
		C /= delta;
		ID = Int(C[0]) * n[1] * n[2] + Int(C[1]) * n[2] + Int(C[2]);
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];
		pmesh->c.push_back(c);
		cLoc[ID].push_back(i);
		blockIndex[i] = ID;
		
		/* mark vertices and facets */
		for(j = 0;j < c.size();j++) {
			Facet& f = mo.f[c[j]];
			(*pfLoc)[c[j]] = 1;
			for(k = 0;k < f.size();k++) {
				(*pvLoc)[f[k]] = 1;	
			}
		}
	}
	/*add vertices & cells*/
	for(ID = 0;ID < total;ID++) {
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];

		count = 0;
		for(i = 0;i < mo.v.size();i++) {
			if((*pvLoc)[i]) {
				pmesh->v.push_back(mo.v[i]);
				(*pvLoc)[i] = count++;
			} else
				(*pvLoc)[i] = Constants::MAX_INT;
		}

		count = 0;
		for(i = 0;i < mo.f.size();i++) {
			if((*pfLoc)[i]) {
				pmesh->f.push_back(mo.f[i]);
				(*pfLoc)[i] = count++;
			} else
				(*pfLoc)[i] = Constants::MAX_INT;
		}
	}
	/*adjust IDs*/
	for(ID = 0;ID < total;ID++) {
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];

		for(i = 0;i < pmesh->f.size();i++) {
			Facet& f = pmesh->f[i];
			for(j = 0;j < f.size();j++)
				f[j] = (*pvLoc)[f[j]];
		}

		for(i = 0;i < pmesh->c.size();i++) {
			Cell& c = pmesh->c[i];
			for(j = 0;j < c.size();j++)
				c[j] = (*pfLoc)[c[j]];
		}
	}
	/*inter mesh faces*/
	IntVector* imesh = new IntVector[total * total];
	Int co,cn;
	for(i = 0;i < mo.f.size();i++) {
		if(gFN[i] != MAX_INT) {
			co = blockIndex[gFO[i]];
			cn = blockIndex[gFN[i]];
			if(co != cn) {
				imesh[co * total + cn].push_back(fLoc[co][i]);
				imesh[cn * total + co].push_back(fLoc[cn][i]);
			}
		}
	}

	/*write meshes to file */
	for(ID = 0;ID < total;ID++) {
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];
		
		/*create directory and switch to it*/
		stringstream path;
		path << mo.name << ID;

		System::mkdir(path.str());
		if(!System::cd(path.str()))    
			return 1;

		/*v,f & c*/
		ofstream of(path.str().c_str());
		of << hex;
		of << pmesh->v << endl;
		of << pmesh->f << endl;
		of << pmesh->c << endl;

		/*bcs*/
		for(Boundaries::iterator it = mo.bdry.begin();it != mo.bdry.end();++it) {
			IntVector b;	
			Int f;
			for(j = 0;j < it->second.size();j++) {
				f = (*pfLoc)[it->second[j]];
				if(f != Constants::MAX_INT)
					b.push_back(f);
			}
			/*write to file*/
			if(b.size()) {
				of << it->first << "  ";
				of << b << endl;
			}
		}

		/*index file*/
		ofstream of2("index");
		of2 << cLoc[ID] << endl;
		
		/*inter mesh boundaries*/
		for(j = 0;j < total;j++) {
			IntVector& f = imesh[ID * total + j];
			if(f.size()) {
				of << "interMesh_" << ID << "_" << j << " ";
				of << f << endl;
				of2 << "interMesh_" << ID << "_" << j << " "
					<< "GHOST" << endl;
			}
		}

		of << dec;
		/*go back*/
		if(!System::cd("..")) 
			return 1;
	}

	/*delete*/
	delete[] meshes;
	delete[] imesh;
	delete[] vLoc;
	delete[] fLoc;
	delete[] cLoc;
	return 0;
}
