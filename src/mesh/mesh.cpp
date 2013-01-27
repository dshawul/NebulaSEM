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
	IntVector&       gFO = gMesh.fo;
	IntVector&       gFN = gMesh.fn;
	Int&             gBCellsStart = gMesh.nc;
	vector<BasicBCondition*> AllBConditions;
	std::vector<interBoundary>& gInterMesh = gMesh.interMesh;
	Vertices         probePoints;

	std::vector<Vector> _fC;
	std::vector<Vector> _cC;
	std::vector<Vector> _fN;
	std::vector<Scalar> _cV;
	std::vector<bool>   _reversed;
}

/*read mesh*/
void Mesh::readMesh() {
	cout << "Reading mesh :" << endl;
	ifstream is(gMeshName.c_str());
	is >> hex;
	is >> gVertices;
	cout << " \t" << gVertices.size() << " vertices" << endl;
	is >> gFacets;
	cout << " \t" << gFacets.size() << " facets" << endl;
	is >> gCells;
	cout << " \t" << gCells.size() << " cells" << endl;
	cout << "Boundaries :" << endl;
	while(Util::nextc(is)) {
		IntVector index;
		string str;
		is >> str;
		cout << " \t" << str << endl;
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
}
/*write mesh*/
void Mesh::MeshObject::write(ostream& os) {
	os << hex;
	os << v;
	os << f;
	os << c;
	for(Boundaries::iterator it = bdry.begin();it != bdry.end();++it)
		os << it->first << " " << it->second << endl;
	os << dec;
}
/*Is face in boundary*/
bool Mesh::faceInBoundary(Int f) {
	for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it) {
		IntVector& gB = it->second;
		forEach(gB,j) {
			if(gB[j] == f)
				return true;
		}
	}
	return false;
}
/*add boundary cells*/
void Mesh::addBoundaryCells() {
	using namespace Constants;
	Int i,index;

	/*neighbor and owner cells of face*/
	gBCellsStart = gCells.size();
	gFO.assign(gFacets.size(),MAX_INT);
	gFN.assign(gFacets.size(),MAX_INT);
	for(i = 0;i < gBCellsStart;i++) {
		forEach(gCells[i],j) {
			index = gCells[i][j];
			if(gFO[index] == MAX_INT) 
				gFO[index] = i;
			else 
				gFN[index] = i;
		}
	}
	/*Flag boundary faces not in gBoundaries for auto deletion*/
	IntVector& gDelete = gBoundaries["delete"];
	forEach(gFN,i) {
		if(gFN[i] == MAX_INT) {
			if(!faceInBoundary(i))
				gDelete.push_back(i);
		}
	}
	/*add boundary cells*/
	for(Boundaries::iterator it = gBoundaries.begin();
		it != gBoundaries.end();++it) {
		IntVector& facets = it->second;
		forEach(facets,j) {
			i = facets[j];
			/*external patch*/
			if(gFN[i] == MAX_INT) {
				Cell c;
				c.push_back(i);
				gCells.push_back(c);
				gFN[i] = gCells.size() - 1;
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
		v = _fC[i] - _cC[gFO[i]];
		if((v & N) < 0) {
			N = -N;
			_reversed[i] = true;
		}
		_fN[i] = N / Scalar(2);
	}
	/* cell volumes */
	for(i = 0;i < gBCellsStart;i++) {
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
	for(i = gBCellsStart;i < gCells.size();i++) {
		_cV[i] = _cV[gFO[gCells[i][0]]];
		_cC[i] = _fC[gCells[i][0]];
	}
}
/* 
* Remove empty boundary
*/
void Mesh::removeBoundary(IntVector& fs) {
	cout << "Removing faces: " << fs.size() << endl;

	Int count;
	IntVector Idf(gFacets.size(),0);
	IntVector Idc(gCells.size(),0);

	/*erase facet reference*/
	forEach(fs,i) {
		Int f = fs[i];
		Cell& co = gCells[gFO[f]];
		forEach(co,j) {
			if(co[j] == f) {
				co.erase(co.begin() + j); 
				break; 
			}
		}
		Cell& cn = gCells[gFN[f]];
		forEach(cn,j) {
			if(cn[j] == f) { 
				cn.erase(cn.begin() + j); 
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
	forEach(gFacets,i) {
		if(gFacets[i].size() == 0) {
			gFacets.erase(gFacets.begin() + i);
			gFO.erase(gFO.begin() + i);
			gFN.erase(gFN.begin() + i);
			_fC.erase(_fC.begin() + i);
			_fN.erase(_fN.begin() + i);
			--i;
		}
	}
	/*updated facet id*/
	count = 0;
	forEach(gCells,i) {
		if(gCells[i].size() != 0) 
			Idc[i] = count++;
		else
			Idc[i] = Constants::MAX_INT;
	}
	/*erase cells*/
	forEach(gCells,i) {
		if(gCells[i].size() == 0) { 
			gCells.erase(gCells.begin() + i); 
			_cC.erase(_cC.begin() + i);
			_cV.erase(_cV.begin() + i);
			--i;
		} else {
			forEach(gCells[i],j) {
				gCells[i][j] = Idf[gCells[i][j]];
			}
		}
	}
	/*facet owner and neighbor*/
	forEach(gFacets,i) {
		gFO[i] = Idc[gFO[i]];
		gFN[i] = Idc[gFN[i]];
	}
	/*patches*/
	for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it) {
		IntVector& gB = it->second;
		forEach(gB,i)
			gB[i] = Idf[gB[i]];
	}

	cout << "Total faces: " << gFacets.size() << endl;
}
/*find nearest cell*/
int Mesh::findNearestCell(const Vector& v) {
	Scalar mindist,dist;
	int bi;
	bi = 0;
	mindist = mag(v - _cC[0]);
	forEach(gCells,i) {
		dist = mag(v - _cC[i]);
		if(dist < mindist) {
			mindist = dist;
			bi = i;
		}
	}
	return bi;
}
