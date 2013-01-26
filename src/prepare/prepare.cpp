#include <sstream>
#include "prepare.h"
#include "system.h"

using namespace std;
using namespace Mesh;

/*duplicate fields*/
template <class T>
void duplicateFields(istream& is,ostream& of) {
	MeshField<T,CELL> f;

	/*internal*/
	f.readInternal(is);

	/*index file*/
	IntVector cLoc;
	ifstream index("index");
	index >> cLoc;

	/*write it out*/
	of << "size "<< sizeof(T) / sizeof(Scalar) << endl;
	of << cLoc.size() << endl;
	of << "{" << endl;
	for(Int j = 0;j < cLoc.size();j++)
		of << f[cLoc[j]] << endl;
	of << "}" << endl;

	/*boundaries*/
	char c;
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
void Prepare::decomposeFields(vector<string>& fields,std::string mName,Int start_index) {
	int size;
	std::string str;

	for(Int ID = start_index;;ID++) {
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
				case 1 :  duplicateFields<Scalar>(is,of); break;
				case 3 :  duplicateFields<Vector>(is,of); break;
				case 6 :  duplicateFields<STensor>(is,of); break;
				case 9 :  duplicateFields<Tensor>(is,of); break;
				}
				/*end*/
			}
		}
		/*go back*/
		System::cd("..");
	}
}
/*decompose in x,y,z direction*/
int Prepare::decomposeXYZ(Mesh::MeshObject& mo,Int* n) {
	
	using Constants::MAX_INT;
	Int i,j,k,ID,count,total = n[0] * n[1] * n[2];
	Vector maxV(Scalar(-10e30)),minV(Scalar(10e30)),delta;

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
	for(j = 0;j < 3;j++) 
		delta[j] /= Scalar(n[j]);

	/*decompose cells*/
	MeshObject *pmesh;
	IntVector *pvLoc,*pfLoc,blockIndex;
	Vector C;
	
	blockIndex.assign(gBCellsStart,0);

	for(i = 0;i < gBCellsStart;i++) {
		Cell& c = mo.c[i];

		/* add cell */
		C = (_cC[i] - minV) / delta;
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
		if(gFN[i] < gBCellsStart) {
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
		ofstream of(mo.name.c_str());
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
					<< "{\n\ttype GHOST\n}" << endl;
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
/*read fields*/
template <class T>
void readFields(istream& is,void* pFields,const IntVector& cLoc) {
	MeshField<T,CELL>& f = *((MeshField<T,CELL>*)pFields);
	Int size;
	char symbol;
	is >> size >> symbol;
	for(Int j = 0;j < size;j++) {
		is >> f[cLoc[j]];
	}
	is >> symbol;
}
/*create fields*/
void createFields(vector<string>& fields,void**& pFields) {
	std::string str;
	Int size;

	/*for each field*/
	pFields = new void*[fields.size()];
	for(Int i = 0;i < fields.size();i++) {
		/*read at time 0*/
		str = fields[i] + "0";
		ifstream is(str.c_str());
		if(!is.fail()) {
			/*fields*/
			is >> str >> size;
			switch(size) {
				case 1 :  pFields[i] = new ScalarCellField(fields[i].c_str(),READWRITE); break;
				case 3 :  pFields[i] = new VectorCellField(fields[i].c_str(),READWRITE); break;
				case 6 :  pFields[i] = new STensorCellField(fields[i].c_str(),READWRITE); break;
				case 9 :  pFields[i] = new TensorCellField(fields[i].c_str(),READWRITE); break;
			}
			/*end*/
		}
	}
}
/*Reverse decomposition*/
int Prepare::merge(Mesh::MeshObject& mo,Int* n,
					 vector<string>& fields,std::string mName,Int start_index) {
    /*create fields*/
	void** pFields;
	createFields(fields,pFields);

	/*indexes*/
	Int total = n[0] * n[1] * n[2];
	IntVector* cLoc = new IntVector[total];
	std::string str;
	Int size;

	for(Int ID = 0;ID < total;ID++) {
		stringstream path;
		path << mName << ID;
		str = path.str() + "/index"; 
		ifstream index(str.c_str());
		index >> cLoc[ID];
	}

	/*for each time step*/
	Int step = start_index;
	Mesh::read_fields(step);
	for(step = start_index + 1;;step++) {
		Int count = 0;
		for(Int ID = 0;ID < total;ID++) {
			stringstream path;
			path << mName << ID;
			for(Int i = 0;i < fields.size();i++) {
				stringstream fpath;
				fpath << fields[i] << step;
				str = path.str() + "/" + fpath.str();

				ifstream is(str.c_str());
				if(is.fail())
					continue;
				count++;
				is >> str >> size;
				switch(size) {
				case 1 :  readFields<Scalar>(is,pFields[i],cLoc[ID]); break;
				case 3 :  readFields<Vector>(is,pFields[i],cLoc[ID]); break;
				case 6 :  readFields<STensor>(is,pFields[i],cLoc[ID]); break;
				case 9 :  readFields<Tensor>(is,pFields[i],cLoc[ID]); break;
				}
			}
		}
		if(count == 0) break;
		Mesh::write_fields(step);
	}

	return 0;
}
/*Convert to VTK format*/
int Prepare::convertVTK(Mesh::MeshObject& mo,vector<string>& fields,Int start_index) {
    /*create fields*/
	void** pFields;
	createFields(fields,pFields);

	/*for each time step*/
	for(Int step = start_index;;step++) {
		Int count = 0;
		for(Int i = 0;i < fields.size();i++) {
			stringstream fpath;
			fpath << fields[i] << step;
			ifstream is(fpath.str().c_str());
			if(is.fail())
				continue;
			count++;
			break;
		}
		if(count == 0) break;
		Mesh::read_fields(step);

		/*write vtk*/
		Vtk::write_vtk(step);
	}

	return 0;
}
/*Probe values at certain locations*/
int Prepare::probe(Mesh::MeshObject& mo,vector<string>& fields,Int start_index) {
	/*probe points*/
	IntVector probes;
	for(Int j = 0;j < Mesh::probePoints.size();j++) {
		Vector v = Mesh::probePoints[j];
		Int index = Mesh::findNearestCell(v);
		probes.push_back(index);
	}
	ofstream of("probes");

    /*create fields*/
	void** pFields;
	createFields(fields,pFields);

	/*for each time step*/
	for(Int step = start_index;;step++) {
		Int count = 0;
		for(Int i = 0;i < fields.size();i++) {
			stringstream fpath;
			fpath << fields[i] << step;
			ifstream is(fpath.str().c_str());
			if(is.fail())
				continue;
			count++;
			break;
		}
		if(count == 0) break;
		Mesh::read_fields(step);

		/*write probes*/
#define WRITE(T) {												\
		T* pf;													\
		for(std::list<T*>::iterator	it = T::fields_.begin();	\
			it != T::fields_.end();++it) {						\
			pf = *it;											\
			of << (*pf)[probes[i]] << " ";						\
		}														\
}
		for(Int i = 0;i < probes.size();i++) {
			of << step << " " << i << " " << 
				Mesh::probePoints[i] << " ";
			WRITE(ScalarCellField);
			WRITE(VectorCellField);
			WRITE(STensorCellField);
			WRITE(TensorCellField);
			of << endl;
		}
#undef WRITE
	}
	/*close*/
	of.close();
	return 0;
}