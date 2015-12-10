#include "mesh.h"

using namespace std;

/*global mesh*/
namespace Mesh {
	MeshObject       gMesh;
	std::string&     gMeshName = gMesh.name;
	Vertices&        gVertices = gMesh.v;
	Facets&          gFacets   = gMesh.f;
	Cells&           gCells    = gMesh.c;
	Boundaries&      gBoundaries = gMesh.bdry;
	IntVector&       gFOC = gMesh.fo;
	IntVector&       gFNC = gMesh.fn;
	Int&             gBCS = gMesh.nc;
	Int&             gBCSI = gMesh.nci;
	vector<BasicBCondition*> AllBConditions;
	std::vector<interBoundary>& gInterMesh = gMesh.interMesh;
	Vertices         probePoints;
	Cells            faceID;

	std::vector<Vector> _fC;
	std::vector<Vector> _cC;
	std::vector<Vector> _fN;
	std::vector<Scalar> _cV;
	std::vector<bool>   _reversed;
}
void Mesh::clear() {
	gMesh.clear();
	Mesh::clearBC();
	probePoints.clear();
}
/*read mesh*/
bool Mesh::readMesh(Int step,bool first, bool coarse) {
	/*open file*/
	stringstream path;
	if(coarse) 
		path << gMeshName;
	else
		path << gMeshName << "_" << step;
	string str = path.str();
	ifstream is(str.c_str());
	if(is.fail()) {
		if(first) {
			str = gMeshName;
			is.open(str.c_str());
		} else
			return false;
	}
	/*read*/
	clear();
	is >> hex;
	is >> gVertices;
	is >> gFacets;
	is >> gCells;
	while(Util::nextc(is)) {
		IntVector index;
		string str;
		is >> str;
		is >> index;

		IntVector& gB = gBoundaries[str];
		gB.insert(gB.begin(),index.begin(),index.end());

		/*internal mesh boundaries*/
		if(str.find("interMesh") != std::string::npos) {
			interBoundary b;
			sscanf(str.c_str(), "interMesh_%x_%x", &b.from,&b.to);
			b.f    = &gBoundaries[str];
			gInterMesh.push_back(b);
		}
	}
	/*start of buffer*/ 
	Int buffer_index = 0;
	forEach(gInterMesh,i) {
		interBoundary& b = gInterMesh[i];
		b.buffer_index = buffer_index;
		buffer_index += b.f->size();
	}
	is >> dec;
	return true;
}
/*write mesh*/
void Mesh::MeshObject::write(ostream& os) {
	os << hex;
	os.precision(12);
	os << v;
	os.precision(6);
	os << f;
	os << c;
	forEachIt(Boundaries,bdry,it)
		os << it->first << " " << it->second << endl;
	os << dec;
}
/*add boundary cells*/
void Mesh::addBoundaryCells() {
	using namespace Constants;
	
	/*neighbor and owner cells of face*/
	gBCS = gCells.size();
	gFOC.assign(gFacets.size(),MAX_INT);
	gFNC.assign(gFacets.size(),MAX_INT);
	forEach(gCells,i) {
		Cell& c = gCells[i];
		forEach(c,j) {
			Int fi = c[j];
			if(gFOC[fi] == MAX_INT) 
				gFOC[fi] = i;
			else 
				gFNC[fi] = i;
		}
	}
	/*Flag boundary faces not in gBoundaries for auto deletion*/
	{
		IntVector faceInB;
		faceInB.assign(gFacets.size(),0);
		forEachIt(Boundaries,gBoundaries,it) {
			IntVector& gB = it->second;	
			forEach(gB,j)
				faceInB[gB[j]] = 1;
		}
	
		IntVector& gDelete = gBoundaries["delete"];
		forEach(gFNC,i) {
			if(gFNC[i] == MAX_INT) {
				if(!faceInB[i])
					gDelete.push_back(i);
			}
		}
	}
	/*reorder cells*/
	{	
		Cells bcs;
		IntVector allbs;
		Int count = 0;
		Int bdry_size = 0;
		
		forEachIt(Boundaries,gBoundaries,it)
			bdry_size += it->second.size();
		bcs.resize(bdry_size);
		allbs.resize(bdry_size);
		forEach(gCells,i) {
			Cell& c = gCells[i];
			forEach(c,j) {
				Int fi = c[j];
				if(gFNC[fi] == MAX_INT) {
					allbs[count] = i;
					bcs[count] = c;
					count++;
					break;
				}
			}
		}
		gBCSI = gBCS - count;
		allbs.resize(count);
		bcs.resize(count);
		erase_indices(gCells,allbs);
		gCells.insert(gCells.end(),bcs.begin(),bcs.end());
	}
	/*add boundary cells*/
	forEachIt(Boundaries,gBoundaries,it) {
		IntVector& gB = it->second;
		forEach(gB,j) {
			Int fi = gB[j];
			if(gFNC[fi] == MAX_INT) {
				Cell c;
				c.push_back(fi);
				gCells.push_back(c);
				gFNC[fi] = gCells.size() - 1;
			}
		}
	}
}
void Mesh::calcGeometry() {
	Int i;
	/*allocate*/
	_fC.assign(gFacets.size(),Vector(0));
	_cC.assign(gCells.size(),Vector(0));
	_fN.assign(gFacets.size(),Vector(0));
	_cV.assign(gCells.size(),Scalar(0));
	_reversed.assign(gFacets.size(),false);

	/* face centre*/
	forEach(gFacets,i) {
		Facet& f = gFacets[i];
		Vector C(0);
		forEach(f,j)
			C += gVertices[f[j]];
		_fC[i] = C / Scalar(f.size());
	}
	/* cell centre */
	forEach(gCells,i) {
		Cell& c = gCells[i];
		Vector C(0);
		forEach(c,j)
			C += _fC[c[j]];
		_cC[i] = C / Scalar(c.size());
	}
	/* face normal */
	Vector v1,v2,v3,v;
	Scalar magN;
	forEach(gFacets,i) {
		Facet& f = gFacets[i];
		Vector N(0),C(0),Ni;
		Scalar Ntot = Scalar(0);
		v1 = _fC[i];
		forEach(f,j) {
			v2 = gVertices[f[j]];
			if(j + 1 == f.size())
				v3 = gVertices[f[0]];
			else
				v3 = gVertices[f[j + 1]];
			Ni = ((v2 - v1) ^ (v3 - v1));
			magN = mag(Ni);
			C += magN * ((v1 + v2 + v3) / 3);
			Ntot += magN;
			N += Ni;
		}
		_fC[i] = C / Ntot;    /*corrected face centre*/
		v = _fC[i] - _cC[gFOC[i]];
		if((v & N) < 0) {
			N = -N;
			_reversed[i] = true;
		}
		_fN[i] = N / Scalar(2);
	}
	/* cell volumes */
	for(i = 0;i < gBCS;i++) {
		Cell& c = gCells[i];
		Scalar V(0),Vi;
		Vector v = _cC[i],C(0);
		forEach(c,j) {
			v = _cC[i] - _fC[c[j]];
			Vi = mag(v & _fN[c[j]]);
			C += Vi * (2 * _fC[c[j]] + _cC[i]) / 3;
			V += Vi;
		}
		_cC[i] = C / V;         /*corrected cell centre */
		_cV[i] = V / Scalar(3);
	}
	/*boundary cell centre and volume*/
	forEachS(gCells,i,gBCS) {
		Int fi = gCells[i][0];
		_cV[i] = _cV[gFOC[fi]];
		_cC[i] = _fC[fi];
	}
	/*facet ids*/
	faceID.clear();
	forEach(gCells,i) {
		IntVector b;
		forEach(gCells[i],j)
			b.push_back(j);
		faceID.push_back(b);
	}
}
/* 
* Remove empty boundary
*/
void Mesh::removeBoundary(IntVector& fs) {
	Int count;
	IntVector Idf(gFacets.size(),0);
	IntVector Idc(gCells.size(),0);
	
	/*erase facet reference*/
	forEach(fs,i) {
		Int f = fs[i];
		Cell& co = gCells[gFOC[f]];
		Cell& coid = faceID[gFOC[f]];
		forEach(co,j) {
			if(co[j] == f) {
				co.erase(co.begin() + j); 
				coid.erase(coid.begin() + j);
				break; 
			}
		}
		Cell& cn = gCells[gFNC[f]];
		Cell& cnid = faceID[gFNC[f]];
		forEach(cn,j) {
			if(cn[j] == f) { 
				cn.erase(cn.begin() + j); 
				cnid.erase(cnid.begin() + j); 
				break; 
			}
		}
	}
	
	/*updated facet id*/
	forEach(fs,i)
		Idf[fs[i]] = Constants::MAX_INT;
	count = 0;
	forEach(gFacets,i) {
		if(Idf[i] != Constants::MAX_INT) 
			Idf[i] = count++;
		else
			gFacets[i].clear();
	}

	/*erase facets*/
	IntVector fzeroIndices;
	forEach(gFacets,i) {
		if(gFacets[i].size() == 0)
			fzeroIndices.push_back(i);
	}
	erase_indices(gFacets,fzeroIndices);
	erase_indices(gFOC,fzeroIndices);
	erase_indices(gFNC,fzeroIndices);
	erase_indices(_fC,fzeroIndices);
	erase_indices(_fN,fzeroIndices);
	/*updated cell id*/
	count = 0;
	forEach(gCells,i) {
		if(gCells[i].size() != 0) 
			Idc[i] = count++;
		else
			Idc[i] = Constants::MAX_INT;
	}
	/*erase cells*/
	IntVector czeroIndices;
	forEach(gCells,i) {
		if(gCells[i].size() == 0)
			czeroIndices.push_back(i);
	}
	erase_indices(gCells,czeroIndices);
	erase_indices(faceID,czeroIndices);
	erase_indices(_cC,czeroIndices);
	erase_indices(_cV,czeroIndices);
	
	/*updated facet id*/
	forEach(gCells,i) {
		forEach(gCells[i],j) {
			gCells[i][j] = Idf[gCells[i][j]];
		}
	}
	/*facet owner and neighbor*/
	forEach(gFacets,i) {
		gFOC[i] = Idc[gFOC[i]];
		gFNC[i] = Idc[gFNC[i]];
	}
	/*patches*/
	forEachIt(Boundaries,gBoundaries,it) {
		IntVector& gB = it->second;
		forEach(gB,i)
			gB[i] = Idf[gB[i]];
	}
}
