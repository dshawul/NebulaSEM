#include <sstream>
#include "prepare.h"
#include "system.h"
#include "metis.h"

using namespace std;
using namespace Mesh;

/**
duplicate fields
*/
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
/**
decompose fields
*/
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

				/*seekg to beginning*/
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
/**
Decompose in x,y,z direction
*/
void Prepare::decomposeXYZ(Int* n,Scalar* nq,IntVector& blockIndex) {
	Int i,j,ID;
	Vector maxV(Scalar(-10e30)),minV(Scalar(10e30)),delta;
	Vector axis(nq[0],nq[1],nq[2]);
	Scalar theta = nq[3];
	Vector C;
	
	/*max and min points*/
	forEach(gVertices,i) {
		C = rotate(gVertices[i],axis,theta);
		for(j = 0;j < 3;j++) {
			if(C[j] > maxV[j]) maxV[j] = C[j];
			if(C[j] < minV[j]) minV[j] = C[j];
		}
	}
    delta = maxV - minV;
	for(j = 0;j < 3;j++) 
		delta[j] /= Scalar(n[j]);
		
	/*assign block indices to cells*/
	for(i = 0;i < gBCellsStart;i++) {
		C = rotate(_cC[i],axis,theta);
		C = (C - minV) / delta;
		ID = Int(C[0]) * n[1] * n[2] + 
			 Int(C[1]) * n[2] + 
			 Int(C[2]);
		blockIndex[i] = ID;
	}
}
/**
Decompose by cell indices
*/
void Prepare::decomposeIndex(Int total,IntVector& blockIndex) {
	for(Int i = 0;i < gBCellsStart;i++)
		blockIndex[i] = (i / (gBCellsStart / total));
}
/**
Decompose using METIS 5.0
*/
void Prepare::decomposeMetis(int total,IntVector& blockIndex) {
	int ncon = 1;
	int edgeCut = 0;
	int ncells = gBCellsStart;
	std::vector<int> xadj,adjncy;
	std::vector<int> options(METIS_NOPTIONS);
	
	/*default options*/
    METIS_SetDefaultOptions(&options[0]);
    
    /*build adjacency*/
    for(Int i = 0;i < gBCellsStart;i++) {
    	Cell& c = gCells[i];
    	xadj.push_back(adjncy.size());
    	forEach(c,j) {
    		Int f = c[j];
			if(i == gFO[f]) {
				if(gFN[f] < gBCellsStart)
					adjncy.push_back(gFN[f]);
			} else {
				if(gFO[f] < gBCellsStart)
					adjncy.push_back(gFO[f]);
			}
		}
    }
    xadj.push_back(adjncy.size());

    /*partition*/
	METIS_PartGraphRecursive (
		&ncells,
		&ncon,
		&xadj[0],
		&adjncy[0],
		NULL,
		NULL,
		NULL,
		&total,
		NULL,
		NULL,
		&options[0],
		&edgeCut,
		(int*)(&blockIndex[0])
    );
}
/**
Decompose
*/
int Prepare::decompose(Int* n,Scalar* nq,int type) {	
	using Constants::MAX_INT;
	Int i,j,ID,count,total = n[0] * n[1] * n[2];

	/*decomposed mesh*/
	MeshObject* meshes = new MeshObject[total];
	IntVector* vLoc = new IntVector[total];
	IntVector* fLoc = new IntVector[total];
	IntVector* cLoc = new IntVector[total];
	for(i = 0;i < total;i++) {
		vLoc[i].assign(gVertices.size(),0);
		fLoc[i].assign(gFacets.size(),0);
	}

	/*decompose cells*/
	MeshObject *pmesh;
	IntVector *pvLoc,*pfLoc,blockIndex;
	blockIndex.assign(gBCellsStart,0);
	
	/*choose*/
	if(type == 0) 
		decomposeXYZ(n,nq,blockIndex);
	else if(type == 1)
		decomposeIndex(total,blockIndex);
	else if(type == 2)
		decomposeMetis(total,blockIndex);
	else; //default -- assigns all to processor 0
	
	/*add cells*/
	for(i = 0;i < gBCellsStart;i++) {
		Cell& c = gCells[i];

		/* add cell */
		ID = blockIndex[i];
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];
		pmesh->c.push_back(c);
		cLoc[ID].push_back(i);
		
		/* mark vertices and facets */
		forEach(c,j) {
			Facet& f = gFacets[c[j]];
			(*pfLoc)[c[j]] = 1;
			forEach(f,k) {
				(*pvLoc)[f[k]] = 1;	
			}
		}
	}
	
	/*add vertices & facets*/
	for(ID = 0;ID < total;ID++) {
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];

		count = 0;
		forEach(gVertices,i) {
			if((*pvLoc)[i]) {
				pmesh->v.push_back(gVertices[i]);
				(*pvLoc)[i] = count++;
			} else
				(*pvLoc)[i] = Constants::MAX_INT;
		}

		count = 0;
		forEach(gFacets,i) {
			if((*pfLoc)[i]) {
				pmesh->f.push_back(gFacets[i]);
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
	forEach(gFacets,i) {
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
		path << gMeshName << ID;

		System::mkdir(path.str());
		if(!System::cd(path.str()))    
			return 1;

		/*v,f & c*/
		ofstream of(gMeshName.c_str());
		of << hex;
		of << pmesh->v << endl;
		of << pmesh->f << endl;
		of << pmesh->c << endl;

		/*bcs*/
		forEachIt(Boundaries,gMesh.bdry,it) {
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
/**
read fields
*/
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
/**
create fields
*/
void createFields(vector<string>& fields,void**& pFields,Int start_index) {
	std::string str;
	Int size;

	/*for each field*/
	pFields = new void*[fields.size()];
	forEach(fields,i) {
		/*read at time 0*/
		stringstream path;
		path << fields[i] << start_index;
		str = path.str(); 

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
/**
open fields
*/
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
/**
Reverse decomposition
*/
int Prepare::merge(Int* n,vector<string>& fields,std::string mName,Int start_index) {
    /*create fields*/
	void** pFields;
	createFields(fields,pFields,start_index);

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
/**
Convert to VTK format
*/
int Prepare::convertVTK(vector<string>& fields,Int start_index) {
    /*create fields*/
	void** pFields;
	createFields(fields,pFields,start_index);

	/*for each time step*/
	for(Int step = start_index;;step++) {
		if(!checkFields(fields,pFields,step))
			break;

		/*write vtk*/
		Vtk::write_vtk(step);
	}

	return 0;
}
/**
Probe values at specified locations
*/
int Prepare::probe(vector<string>& fields,Int start_index) {
	/*probe points*/
	IntVector probes;
	getProbeFaces(probes);
	ofstream of("probes");

    /*create fields*/
	void** pFields;
	createFields(fields,pFields,start_index);

	/*for each time step*/
	for(Int step = start_index;;step++) {
		if(!checkFields(fields,pFields,step))
			break;

		/*Interpolate*/
		forEachField(interpolateVertexAll());
		
		/*write probes*/
#define ADD(v,value,weight) {										\
		dist = magSq((v) - probeP);									\
		dist = weight / (dist + 1.0f);								\
		sum += (value) * dist;										\
		sumd += dist;												\
}
#define SUM(X) {													\
		Cell& c = gCells[X];										\
		forEach(c,m) {												\
			Facet& f = gFacets[c[m]];								\
			forEach(f,j) {											\
				ADD(gVertices[f[j]],(*it)[f[j]],1.0);				\
			}														\
		}															\
}
#define WRITE(T) {									                \
	    std::list<MeshField<T,CELL>*>::iterator it1 =				\
			MeshField<T,CELL>::fields_.begin();						\
	    for(MeshField<T,CELL>::vertexFieldsType::iterator it =      \
		    (MeshField<T,CELL>::vf_fields_)->begin(); it !=         \
			(MeshField<T,CELL>::vf_fields_)->end(); ++it,++it1) {   \
			T sum(0.0);												\
			Scalar sumd(0.0);										\
			ADD(cC[c1],(*(*it1))[c1],2.0);							\
			ADD(cC[c2],(*(*it1))[c2],2.0);							\
			SUM(sc); 									            \
		    of <<  (sum/sumd) << " ";								\
		}												            \
}
		forEach(probes,i) {
			Int fi = probes[i];
			Int c1 = gFO[fi];
			Int c2 = gFN[fi];
			Vector probeP = probePoints[i];
			Scalar dir = ((fC[fi] - probeP) & fN[fi]),dist;
			Int sc;
			if(dir >= 0) sc = c1;
			else sc = c2;

			of << step << " " << i << " " << probePoints[i] << " ";

			WRITE(Scalar);
			WRITE(Vector);
			WRITE(STensor);
			WRITE(Tensor);

			of << endl;
		}
#undef WRITE
#undef SUM
#undef ADD
	}

	return 0;
}
/**
Refine mesh
*/
void Prepare::refineMesh(const IntVector& rCells, const Int refine_type,
						 const Int refine_shape, const Vector& refine_dir) {

	/*Remove boundary cells*/
	gCells.erase(gCells.begin() + gBCellsStart,gCells.end());
	
	/*build refine flags array*/
	IntVector refineC,refineF,rFacets;
	refineC.assign(gBCellsStart,0);
	refineF.assign(gFacets.size(),0);
	forEach(rCells,i) {
		Int ci = rCells[i];
		refineC[ci] = 1;
		Cell& c = gCells[ci];
		forEach(c,j)
			refineF[c[j]] = 1;
	}
	
	/*choose faces to refine*/
	bool refine3D = equal(Scalar(0.0),mag(refine_dir));
	forEach(refineF,i) {
		if(refineF[i]) {
			Vector eu = unit(fN[i]);
			if(refine3D || 
				equal(refine_dir,eu) || 
				equal(refine_dir,-eu))
				rFacets.push_back(i);
			else
				refineF[i] = 0;
		}
	}
	
	/*refine facets*/
	Int nj;
	IntVector startF;
	startF.assign(gFacets.size(),0);
	Int iBegin = gVertices.size();
	forEach(rFacets,i) {
		Int fi = rFacets[i];
		Facet f = gFacets[fi];
		startF[fi] = gFacets.size();
		Vector C = fC[fi];
		gVertices.push_back(C);
		Int fci = gVertices.size() - 1;
		
		/*quadrilaterals*/
		if(refine_shape == 0) {
			/*add vertices*/
			IntVector midpts;
			Vector v1,v2;
			forEach(f,j) {
				v1 = gVertices[f[j]];
				nj = j + 1;
				if(nj == f.size())
					nj = 0;
				v2 = gVertices[f[nj]];
				Vector Ce = (v1 + v2) / 2.0;
			
				//duplicate
				Int k = iBegin;
				for(; k < gVertices.size();k++) {
					if(equal(Ce,gVertices[k])) {
						midpts.push_back(k);
						break;
					}
				}
				if(k == gVertices.size()) {
					gVertices.push_back(Ce);
					midpts.push_back(k);
				}
			}
		
			/*add facets*/
			Int k1,k2;
			forEach(f,j) {
				k1 = midpts[j];
				if(j == 0)
					nj = f.size() - 1;
				else
					nj = j - 1;
				k2 = midpts[nj];
				/*quad face*/
				Facet fn;
				fn.push_back(f[j]);
				fn.push_back(k2);
				fn.push_back(fci);
				fn.push_back(k1);
				gFacets.push_back(fn);
			}
		/*triangles*/
		} else {
			/*add facets*/
			Int k1,k2;
			forEach(f,j) {
				k1 = f[j];
				nj = j + 1;
				if(nj == f.size())
					nj = 0;
				k2 = f[nj];
				/*quad face*/
				Facet fn;
				fn.push_back(k1);
				fn.push_back(fci);
				fn.push_back(k2);
				gFacets.push_back(fn);
			}
		}
	}
	/*add cells*/
	if(refine_type == 0) {
		forEach(rCells,i) {
			Int ci = rCells[i];
			Cell c = gCells[ci];
			Vector C = cC[ci];
			gVertices.push_back(C);
			Int cci = gVertices.size() - 1;
			Int iBegin = gFacets.size();
			forEach(c,j) {
				Int fi = c[j];
				
				/*cell is refined but face is not?*/
				IntVector list;
				if(!refineF[fi]) {
					list.push_back(fi);
				} else {
					forEach(gFacets[fi],k)
						list.push_back(startF[fi] + k);
				}

				/*refine cells*/
				forEach(list,k) {
					Int fni = list[k];
					Facet f = gFacets[fni];
				
					Cell cn;
					cn.push_back(fni);

					Int v1i,v2i,nj;
					forEach(f,l) {
						v1i = f[l];
						nj = l + 1;
						if(nj == f.size())
							nj = 0;
						v2i = f[nj];
						//triangular face
						Facet fn;
						fn.push_back(v1i);
						fn.push_back(v2i);
						fn.push_back(cci);
						//duplicate
						Int k = iBegin;
						for(; k < gFacets.size();k++) {
							if(equal(fn,gFacets[k])) {
								cn.push_back(k);
								break;
							}
						}
						if(k == gFacets.size()) {
							gFacets.push_back(fn);
							cn.push_back(k);
						}
					}
					gCells.push_back(cn);
				}
			}
		}
	}
	/*adjust facet indexes in cells and boundaries*/
	refineF.resize(gFacets.size(),0);
	Int l = 0;
	forEach(refineF,i) {
		if(!refineF[i]) 
			refineF[i] = l++;
		else
			refineF[i] = Constants::MAX_INT;
	}

	/*erase old facets and add the new ones*/
#define ERASEADD() {\
	IntVector newF,eraseF;\
	forEach(c,j) {\
		Int fi = c[j];\
		if(refineF[fi] == Constants::MAX_INT) {\
			eraseF.push_back(j);\
			forEach(gFacets[fi],k) {\
				Int fni = startF[fi] + k;\
				fni = refineF[fni];\
				newF.push_back(fni);\
			}\
		}\
		c[j] = refineF[fi];\
	}\
	erase_indices(c,eraseF);\
	c.insert(c.end(),newF.begin(),newF.end());\
}

	forEach(gCells,i) {
		Cell& c = gCells[i];
		ERASEADD();
	}
	forEachIt(Boundaries,gBoundaries,it) {
		IntVector& c = it->second;
		ERASEADD();
	}
#undef ERASEADD
	erase_indices(gFacets,rFacets);
	
	/*erase cells*/
	if(refine_type == 0)
		erase_indices(gCells,rCells);
	
	/* Break edge of faces that are not set for refinement 
	 * but should be due to neighboring refined faces*/
	if(refine_shape == 0) {
		forEach(gFacets,i) {
			Facet& f = gFacets[i];
			IntVector addf;
			Vector v1,v2;
			addf.assign(f.size(),Constants::MAX_INT);
			forEach(f,j) {
				v1 = gVertices[f[j]];
				nj = j + 1;
				if(nj == f.size())
					nj = 0;
				v2 = gVertices[f[nj]];
				Vector Ce = (v1 + v2) / 2.0;
		
				//duplicate
				Int k = iBegin;
				for(; k < gVertices.size();k++) {
					if(equal(Ce,gVertices[k])) {
						addf[j] = k;
						break;
					}
				}
			}
			if(addf.size()) {
				Facet nf;
				forEach(f,j) {
					nf.push_back(f[j]);
					if(addf[j] != Constants::MAX_INT)
						nf.push_back(addf[j]);
				}
				f = nf;
			}
		}
	}
	
	/*write*/
	ofstream os(gMeshName.c_str());
	gMesh.write(os);
}