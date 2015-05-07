#include <sstream>
#include "prepare.h"
#include "system.h"
#include "metis.h"

using namespace std;
using namespace Mesh;

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
void createFields(vector<string>& fields,IntVector& field_sizes,void**& pFields,Int step) {
	std::string str;
	Int size;

	/*for each field*/
	pFields = new void*[fields.size()];
	field_sizes.assign(fields.size(),0);
	forEach(fields,i) {
		/*read at time 0*/
		stringstream path;
		path << fields[i] << step;
		str = path.str(); 

		ifstream is(str.c_str());
		if(!is.fail()) {
			/*fields*/
			is >> str >> size;
			field_sizes[i] = size;
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
destroy fields
*/
void destroyFields(void**& pFields,IntVector& field_sizes) {
#define destroy(T) {\
	T x = (T)pFields[i];\
	x->deallocate(false);\
	delete x;\
}
	forEach(field_sizes,i) {
		switch(field_sizes[i]) {
			case 1 :  destroy(ScalarCellField*); break;
			case 3 :  destroy(VectorCellField*); break;
			case 6 :  destroy(STensorCellField*); break;
			case 9 :  destroy(TensorCellField*); break;
		}
	}
#undef destroy
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
		read_fields(step);
	return count;
}
/**
duplicate fields
*/
template <class T>
void duplicateFields(ostream& of,void* pFields,IntVector& cLoc) {
	MeshField<T,CELL>& f = *((MeshField<T,CELL>*)pFields);
	
	/*write it out*/
	of << "size "<< sizeof(T) / sizeof(Scalar) << endl;
	of << cLoc.size() << endl;
	of << "{" << endl;
	forEach(cLoc,j)
		of << f[cLoc[j]] << endl;
	of << "}" << endl;

	/*boundaries*/
	BasicBCondition* bbc;
	BCondition<T>* bc;
	forEach(AllBConditions,i) {
		bbc = AllBConditions[i];
		if(bbc->fIndex == f.fIndex) {
			bc = static_cast<BCondition<T>*> (bbc);
			of << *bc << std::endl;
		}
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
	for(i = 0;i < gBCS;i++) {
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
	for(Int i = 0;i < gBCS;i++)
		blockIndex[i] = (i / (gBCS / total));
}
/**
Decompose using METIS 5.0
*/
void Prepare::decomposeMetis(int total,IntVector& blockIndex) {
	int ncon = 1;
	int edgeCut = 0;
	int ncells = gBCS;
	std::vector<int> xadj,adjncy;
	std::vector<int> options(METIS_NOPTIONS);
	
	/*default options*/
    METIS_SetDefaultOptions(&options[0]);
    
    /*build adjacency*/
    for(Int i = 0;i < gBCS;i++) {
    	Cell& c = gCells[i];
    	xadj.push_back(adjncy.size());
    	forEach(c,j) {
    		Int f = c[j];
			if(i == gFOC[f]) {
				if(gFNC[f] < gBCS)
					adjncy.push_back(gFNC[f]);
			} else {
				if(gFOC[f] < gBCS)
					adjncy.push_back(gFOC[f]);
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
int Prepare::decompose(vector<string>& fields,Int* n,Scalar* nq,int type, Int step) {	
	using Constants::MAX_INT;
	Int i,j,ID,count,total = n[0] * n[1] * n[2];
	
	/*Read mesh*/
	LoadMesh(step,true,false);

	/*Read fields*/
	void** pFields = 0;
	IntVector field_sizes;
	destroyFields(pFields,field_sizes);
	createFields(fields,field_sizes,pFields,step);
	if(!checkFields(fields,pFields,step))
		return 1;
		
	/**********************
	 * decompose mesh
	 **********************/
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
	blockIndex.assign(gBCS,0);
	
	/*choose*/
	if(type == 0) 
		decomposeXYZ(n,nq,blockIndex);
	else if(type == 1)
		decomposeIndex(total,blockIndex);
	else if(type == 2)
		decomposeMetis(total,blockIndex);
	else; //default -- assigns all to processor 0
	
	/*add cells*/
	for(i = 0;i < gBCS;i++) {
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
		if(gFNC[i] < gBCS) {
			co = blockIndex[gFOC[i]];
			cn = blockIndex[gFNC[i]];
			if(co != cn) {
				imesh[co * total + cn].push_back(fLoc[co][i]);
				imesh[cn * total + co].push_back(fLoc[cn][i]);
			}
		}
	}
	/***************************
	 * Expand cLoc
	 ***************************/
	for(ID = 0;ID < total;ID++) {
		const Int block = DG::NP;
		IntVector& cF = cLoc[ID];
		cF.resize(cF.size() * block);
		for(int i = cF.size() - 1;i >= 0;i -= block) {
			Int ii = i / block;
			Int C = cF[ii] * block;
			for(Int j = 0; j < block;j++)
				cF[i - j] = C + block - 1 - j; 
		}
	}
	/***************************
	 * write mesh/index/fields
	 ***************************/
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
		
		/*inter mesh boundaries*/
		for(j = 0;j < total;j++) {
			IntVector& f = imesh[ID * total + j];
			if(f.size()) {
				of << "interMesh_" << ID << "_" << j << " ";
				of << f << endl;
			}
		}
		
		of << dec;
		
		/*index file*/
		stringstream path1;
		path1 << "index_" << step;
		ofstream of2(path1.str().c_str());
		of2 << cLoc[ID] << endl;
		
		/*fields*/
		forEach(fields,i) {
			stringstream path;
			path << fields[i] << step;
			string str = path.str();
			ofstream of3(str.c_str());

			/*fields*/
			switch(field_sizes[i]) {
			case 1 :  duplicateFields<Scalar>(of3,pFields[i],cLoc[ID]); break;
			case 3 :  duplicateFields<Vector>(of3,pFields[i],cLoc[ID]); break;
			case 6 :  duplicateFields<STensor>(of3,pFields[i],cLoc[ID]); break;
			case 9 :  duplicateFields<Tensor>(of3,pFields[i],cLoc[ID]); break;
			}
			
			/*inter mesh boundaries*/
			for(j = 0;j < total;j++) {
				IntVector& f = imesh[ID * total + j];
				if(f.size()) {
					of3 << "interMesh_" << ID << "_" << j << " "
						<< "{\n\ttype GHOST\n}" << endl;
				}
			}
		}
		
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
Reverse decomposition
*/
int Prepare::merge(Int* n,vector<string>& fields,Int start_index) {
	/*indexes*/
	Int total = n[0] * n[1] * n[2];
	IntVector* cLoc = new IntVector[total];
	string str;
	Int size;

	/*merge at each time step*/
	void** pFields = 0;
	IntVector field_sizes;
	for(Int step = start_index;;step++) {
	
		/*load mesh and fields*/
		if(LoadMesh(step,(step == start_index),true)) {
			destroyFields(pFields,field_sizes);
			createFields(fields,field_sizes,pFields,step);
			checkFields(fields,pFields,step);
			for(Int ID = 0;ID < total;ID++) {
				stringstream path;
				path << gMeshName << ID << "/index_" << step;
				ifstream index(path.str().c_str());
				cLoc[ID].clear();
				index >> cLoc[ID];
			}
		}
		if(step == start_index) continue;
		
		/*read and merge fields*/
		Int count = 0;
		for(Int ID = 0;ID < total;ID++) {
			stringstream path;
			path << gMeshName << ID;
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
		write_fields(step);
	}

	return 0;
}
/**
Convert to VTK format
*/
int Prepare::convertVTK(vector<string>& fields,Int start_index) {
	void** pFields = 0;
	IntVector field_sizes;
	for(Int step = start_index;;step++) {
		if(LoadMesh(step,(step == start_index),true)) {
			destroyFields(pFields,field_sizes);
			createFields(fields,field_sizes,pFields,step);
		}
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
	ofstream of("probes");

    /*probe at each time step*/
	void** pFields = 0;
	IntVector field_sizes;
	for(Int step = start_index;;step++) {
		if(LoadMesh(step,(step == start_index),true)) {
			destroyFields(pFields,field_sizes);
			createFields(fields,field_sizes,pFields,step);
			probes.clear();
			getProbeFaces(probes);
		}
		if(!checkFields(fields,pFields,step))
			break;

		/*Interpolate*/
		forEachCellField(interpolateVertexAll());
		
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
refine fields
*/
template <class T>
void refineField(MeshField<T,CELL>& fo,Int step,IntVector& cellMap) {
	MeshField<T,CELL> f(fo.fName.c_str(),fo.access,false);
	forEach(f,i)
		f[i] = fo[cellMap[i]];
	f.write(step);
}
/**
Refine mesh
*/
void Prepare::refineMesh(vector<string>& fields,const RefineParams& rparams, Int step) {
			
	/*Read mesh*/
	LoadMesh(step,true,false);

	/*Read fields*/
	void** pFields = 0;
	IntVector field_sizes;
	destroyFields(pFields,field_sizes);
	createFields(fields,field_sizes,pFields,step);
	if(!checkFields(fields,pFields,step))
		return;
		
	/*find quantity of interest*/
	ScalarCellField qoi;
	vector<string>::iterator it = 
		find(fields.begin(), fields.end(), rparams.field);
	if(it != fields.end()) {
		Int i = it - fields.begin();
		switch(field_sizes[i]) {
		case 1 :  qoi = mag(*((ScalarCellField*)pFields[i])); break;
		case 3 :  qoi = mag(*((VectorCellField*)pFields[i])); break;
		}
	}
	qoi = mag(gradi(qoi));
	Scalar maxq = 0,minq = 10e30;
	for(Int i = 0;i < gBCS;i++) {
		if(qoi[i] > maxq) maxq = qoi[i];
		if(qoi[i] < minq) minq = qoi[i];
	}
	maxq = max(Constants::MachineEpsilon,maxq);
	qoi /= maxq;
	maxq = 1;
	minq = minq / maxq;

	/*get cells to refine*/
	gCells.erase(gCells.begin() + gBCS,gCells.end());
	IntVector rCells;
	for(Int i = 0;i < gBCS;i++) {
		if(qoi[i] >= rparams.field_max)
			rCells.push_back(i);
	}
	/*refine mesh and fields*/
	IntVector cellMap;
	refineMesh(rCells,rparams,cellMap);
	forEach(fields,i) {
		switch(field_sizes[i]) {
		case 1 :  refineField<Scalar>(*((ScalarCellField*)pFields[i]),step,cellMap); break;
		case 3 :  refineField<Vector>(*((VectorCellField*)pFields[i]),step,cellMap); break;
		case 6 :  refineField<STensor>(*((STensorCellField*)pFields[i]),step,cellMap); break;
		case 9 :  refineField<Tensor>(*((TensorCellField*)pFields[i]),step,cellMap); break;
		}
	}
	
	/*Write mesh*/
	stringstream path;
	path << gMeshName << "_" << step;
	ofstream os(path.str().c_str());
	gMesh.write(os);
}

void Prepare::refineMesh(IntVector& rCells,const RefineParams& rparams,IntVector& cellMap) {
			
	/*build refine flags array*/
	IntVector refineC,refineF,
			  rFacets,rsFacets;
	cellMap.assign(gBCS,0);
	forEach(cellMap,i)
		cellMap[i] = i;
	refineC.assign(gBCS,0);
	refineF.assign(gFacets.size(),0);
	forEach(rCells,i) {
		Int ci = rCells[i];
		refineC[ci] = 1;
		Cell& c = gCells[ci];
		forEach(c,j)
			refineF[c[j]] = 1;
	}
	
	/*choose faces to refine*/
	bool refine3D = equal(Scalar(0.0),mag(rparams.dir));
	forEach(refineF,i) {
		if(refineF[i]) {
			Vector eu = unit(fN[i]);
			if(refine3D || 
				equal(rparams.dir,eu) || 
				equal(rparams.dir,-eu))
				rFacets.push_back(i);
			else
				refineF[i] = 0;
		}
	}
	
	/***************************
	 * Refine facets
	 ***************************/
	Int nj;
	IntVector startF;
	startF.assign(gFacets.size(),0);
	Int ivBegin = gVertices.size();
	forEach(rFacets,i) {
		Int fi = rFacets[i];
		Facet f = gFacets[fi];
		startF[fi] = gFacets.size();
		Vector C = fC[fi];
		gVertices.push_back(C);
		Int fci = gVertices.size() - 1;
		
		/*quadrilaterals*/
		if(rparams.shape == 0) {
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
				Int k = ivBegin;
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
	/*************************
	 * Refine cells
	 *************************/
	forEach(rCells,i) {
		Int ci = rCells[i];
		Cell c = gCells[ci];
		Vector C = cC[ci];
		gVertices.push_back(C);
		Int cci = gVertices.size() - 1;
		Int ifBegin = gFacets.size();
		Cells newc;
		
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
					Int k = ifBegin;
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
				newc.push_back(cn);
			}
		}
		/*merge cells*/
		if(rparams.shape == 0) {
			/*cells on same face*/
			IntVector ownerf;
			forEach(c,j) {
				Int fi = c[j];
				forEach(gFacets[fi],k)
					ownerf.push_back(fi);
			}
			/*find cells to merge with a shared face*/
			Cells mergec;
			forEach(newc,j) {
				Int o1 = ownerf[j];
				Cell& c1 = newc[j];
				IntVector mg;
				forEachS(newc,k,j+1) {
					Int o2 = ownerf[k];
					Cell& c2 = newc[k];
					if(o1 == o2) continue;
					forEach(c1,m) {
						Int f1 = c1[m];
						forEach(c2,n) {
							Int f2 = c2[n];
							if(f1 == f2) {
								mg.push_back(k);
								c1.erase(c1.begin() + m);
								c2.erase(c2.begin() + n);
								rsFacets.push_back(f1);
								goto END;
							}
						}
					}
					END:;					
				}
				mergec.push_back(mg);
			}
			/*merge cells*/
			IntVector erasei;
			erasei.assign(newc.size(),0);
			forEachRev(mergec,j) {
				IntVector& cm = mergec[j];
				forEach(cm,k)
					erasei[cm[k]] = 1;
				if(cm.size()) {
					Cell& c1 = newc[j];
					Cell& c2 = newc[cm[0]];
					c1.insert(c1.end(),c2.begin(),c2.end());
				}
			}
			/*erase cells*/
			IntVector erasec;
			forEach(erasei,j) {
				if(erasei[j])
					erasec.push_back(j);
			}
			erase_indices(newc,erasec);
			/*two cells should share only one face*/
			forEach(newc,j) {
				Cell& c1 = newc[j];
				forEachS(newc,k,j+1) {
					Cell& c2 = newc[k];
					//find shared faces
					IntVector shared1,shared2;
					forEach(c1,m) {
						Int f1 = c1[m];
						forEach(c2,n) {
							Int f2 = c2[n];
							if(f1 == f2) {
								shared1.push_back(m);
								shared2.push_back(n);
							}
						}
					}
					//two faces shared between two cells
					if(shared1.size() == 2) {
						Int fi1 = c1[shared1[0]];
						Int fi2 = c1[shared1[1]];
						Facet& f1 = gFacets[fi1];
						Facet& f2 = gFacets[fi2];
						erase_indices(c1,shared1);
						erase_indices(c2,shared2);
						rsFacets.push_back(fi1);
						rsFacets.push_back(fi2);
						//add new face by merging two faces
						Facet f;
						if(f1[0] == f2[0]) {
							f.push_back(f1[1]);
							f.push_back(f1[0]);
							f.push_back(f2[1]);
							f.push_back(f2[2]);
						} else if(f1[0] == f2[1]) {
							f.push_back(f1[1]);
							f.push_back(f1[0]);
							f.push_back(f2[0]);
							f.push_back(f2[2]);
						} else if(f1[1] == f2[1]) {
							f.push_back(f1[0]);
							f.push_back(f1[1]);
							f.push_back(f2[0]);
							f.push_back(f2[2]);
						} else if(f1[1] == f2[0]) {
							f.push_back(f1[0]);
							f.push_back(f1[1]);
							f.push_back(f2[1]);
							f.push_back(f2[2]);
						}
						gFacets.push_back(f);
						Int fi = gFacets.size() - 1;
						c1.push_back(fi);
						c2.push_back(fi);
					}
				}
			}
		}
		/*add cells*/
		gCells.insert(gCells.end(),newc.begin(),newc.end());
		forEach(newc,j)
			cellMap.push_back(ci);
	}
	/*************************************
	 *  Remove refined facets and cells
	 *************************************/
	/*add additional facets to remove*/
	refineF.resize(gFacets.size(),0);
	forEach(rsFacets,i) {
		Int fi = rsFacets[i];
		if(!refineF[fi]) {
			refineF[fi] = Constants::MAX_INT - 1;
			rFacets.push_back(fi);
		}
	}
	/*adjust facet indexes in cells and boundaries*/
	Int count = 0;
	forEach(refineF,i) {
		if(refineF[i] == 0) refineF[i] = count++;
		else if(refineF[i] == Constants::MAX_INT - 1);
		else refineF[i] = Constants::MAX_INT;
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

	/*erase facets*/
	sort(rFacets.begin(),rFacets.end());
	erase_indices(gFacets,rFacets);
	
	/*erase cells*/
	sort(rCells.begin(),rCells.end());
	erase_indices(gCells,rCells);
	erase_indices(cellMap,rCells);
	gBCS = gCells.size();
	gBCSfield = gBCS * DG::NP;
	
	/*******************************************************
	 * Break edge of faces that are not set for refinement 
	 * but should be due to neighboring refined faces 
	 *******************************************************/
	if(rparams.shape == 0) {
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
				Int k = ivBegin;
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
}
