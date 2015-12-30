#include "prepare.h"
#include "system.h"
#include "metis.h"

using namespace std;
using namespace Mesh;

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
		C = rotate(gCC[i],axis,theta);
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
	
	/*no decomposition*/
	if(total == 1)
		return 1;
	
	/*Read mesh*/
	LoadMesh(step,true,false);

	/*Read fields*/
	createFields(fields,step);
	if(!readFields(fields,step))
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
		pmesh->mCells.push_back(c);
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
				pmesh->mVertices.push_back(gVertices[i]);
				(*pvLoc)[i] = count++;
			} else
				(*pvLoc)[i] = Constants::MAX_INT;
		}

		count = 0;
		forEach(gFacets,i) {
			if((*pfLoc)[i]) {
				pmesh->mFacets.push_back(gFacets[i]);
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

		forEach(pmesh->mFacets,i) {
			Facet& f = pmesh->mFacets[i];
			forEach(f,j)
				f[j] = (*pvLoc)[f[j]];
		}

		forEach(pmesh->mCells,i) {
			Cell& c = pmesh->mCells[i];
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
	char wdir[512]; 
	System::pwd(wdir,512);
	
	for(ID = 0;ID < total;ID++) {
		pmesh = &meshes[ID];
		pvLoc = &vLoc[ID];
		pfLoc = &fLoc[ID];
		
		/*create directory and switch to it*/
		stringstream path;
		path << gMeshName << ID;

		System::cd(wdir);
		System::mkdir(path.str());
		if(!System::cd(path.str()))    
			return 1;

		/*v,f & c*/
		ofstream of(gMeshName.c_str());
		of << hex;
		of << pmesh->mVertices << endl;
		of << pmesh->mFacets << endl;
		of << pmesh->mCells << endl;

		/*bcs*/
		forEachIt(Boundaries,gMesh.mBoundaries,it) {
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
			BaseField* pf = BaseField::findField(fields[i]);
			pf->writeInternal(of3,&cLoc[ID]);
			pf->writeBoundary(of3);

			/*inter mesh boundaries*/
			for(j = 0;j < total;j++) {
				IntVector& f = imesh[ID * total + j];
				if(f.size()) {
					of3 << "interMesh_" << ID << "_" << j << " "
						<< "{\n\ttype GHOST\n}" << endl;
				}
			}
		}
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

	/*merge at each time step*/
	for(Int step = start_index;;step++) {
	
		/*load mesh and fields*/
		if(LoadMesh(step,(step == start_index),true)) {
			createFields(fields,step);
			readFields(fields,step);
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
				BaseField* pf = BaseField::findField(fields[i]);
				pf->readInternal(is,&cLoc[ID]);
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
	for(Int step = start_index;;step++) {
		if(LoadMesh(step,(step == start_index),true))
			createFields(fields,step);
		if(!readFields(fields,step))
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
	for(Int step = start_index;;step++) {
		if(LoadMesh(step,(step == start_index),true)) {
			createFields(fields,step);
			probes.clear();
			getProbeFaces(probes);
		}
		if(!readFields(fields,step))
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
			Int c1 = FO[fi];
			Int c2 = FN[fi];
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


