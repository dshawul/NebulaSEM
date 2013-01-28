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
	forEach(cLoc,j)
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
		forEach(fields,i) {
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
int Prepare::decomposeXYZ(Mesh::MeshObject& mo,Int* n,Scalar* nq) {
	
	using Constants::MAX_INT;
	Int i,j,ID,count,total = n[0] * n[1] * n[2];
	Vector maxV(Scalar(-10e30)),minV(Scalar(10e30)),delta;
	Vector axis(nq[0],nq[1],nq[2]);
	Scalar theta = nq[3];
	Vector C;

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
	forEach(mo.v,i) {
		C = rotate(mo.v[i],axis,theta);
		for(j = 0;j < 3;j++) {
			if(C[j] > maxV[j]) maxV[j] = C[j];
			if(C[j] < minV[j]) minV[j] = C[j];
		}
	}
    delta = maxV - minV;
	for(j = 0;j < 3;j++) 
		delta[j] /= Scalar(n[j]);

	/*decompose cells*/
	MeshObject *pmesh;
	IntVector *pvLoc,*pfLoc,blockIndex;
	
	blockIndex.assign(gBCellsStart,0);

	for(i = 0;i < gBCellsStart;i++) {
		Cell& c = mo.c[i];

		/* add cell */
		C = rotate(_cC[i],axis,theta);
		C = (C - minV) / delta;
		ID = Int(C[0]) * n[1] * n[2] + 
			 Int(C[1]) * n[2] + 
			 Int(C[2]);
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];
		pmesh->c.push_back(c);
		cLoc[ID].push_back(i);
		blockIndex[i] = ID;
		
		/* mark vertices and facets */
		forEach(c,j) {
			Facet& f = mo.f[c[j]];
			(*pfLoc)[c[j]] = 1;
			forEach(f,k) {
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
		forEach(mo.v,i) {
			if((*pvLoc)[i]) {
				pmesh->v.push_back(mo.v[i]);
				(*pvLoc)[i] = count++;
			} else
				(*pvLoc)[i] = Constants::MAX_INT;
		}

		count = 0;
		forEach(mo.f,i) {
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

		forEach(pmesh->f,i) {
			Facet& f = pmesh->f[i];
			forEach(f,j)
				f[j] = (*pvLoc)[f[j]];
		}

		forEach(pmesh->c,i) {
			Cell& c = pmesh->c[i];
			forEach(c,j)
				c[j] = (*pfLoc)[c[j]];
		}
	}
	/*inter mesh faces*/
	IntVector* imesh = new IntVector[total * total];
	Int co,cn;
	forEach(mo.f,i) {
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
		forEachIt(Boundaries,mo.bdry,it) {
			IntVector b;	
			Int f;
			forEach(it->second,j) {
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
	forEach(fields,i) {
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
/*open fields*/
Int checkFields(vector<string>& fields,void**& pFields,Int step) {
	Int count = 0;
	forEach(fields,i) {
		stringstream fpath;
		fpath << fields[i] << step;
		ifstream is(fpath.str().c_str());
		if(is.fail())
			continue;
		count++;
		break;
	}
	if(count)
		Mesh::read_fields(step);
	return count;
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
			forEach(fields,i) {
				stringstream fpath;
				fpath << fields[i] << step;
				str = path.str() + "/" + fpath.str();
				ifstream is(str.c_str());
				if(is.fail())
					continue;
				count++;
				/*read*/
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
		if(!checkFields(fields,pFields,step))
			break;

		/*write vtk*/
		Vtk::write_vtk(step);
	}

	return 0;
}
/*Probe values at certain locations*/
int Prepare::probe(Mesh::MeshObject& mo,vector<string>& fields,Int start_index) {
	/*probe points*/
	IntVector probes;
	getProbeCells(probes);
	ofstream of("probes");

    /*create fields*/
	void** pFields;
	createFields(fields,pFields);

	/*for each time step*/
	for(Int step = start_index;;step++) {
		if(!checkFields(fields,pFields,step))
			break;

		/*write probes*/
#define WRITE(T) {												\
		forEachIt(std::list<T*>,T::fields_,it)  				\
			of << (*(*it))[probes[i]] << " ";					\
}
		forEach(probes,i) {
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

	return 0;
}