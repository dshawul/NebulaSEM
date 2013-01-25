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
	for(Int i = 0;i < gInterMesh.size();i++) {
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
		for(Int j = 0;j < gB.size();j++) {
			if(gB[j] == f) {
				return true;
			}
		}
	}
	return false;
}
/*add boundary cells*/
void Mesh::addBoundaryCells() {
	using namespace Constants;
	Int i,j,index;

	/*neighbor and owner cells of face*/
	gBCellsStart = gCells.size();
	gFO.assign(gFacets.size(),MAX_INT);
	gFN.assign(gFacets.size(),MAX_INT);
	for(i = 0;i < gBCellsStart;i++) {
		for(j = 0;j < gCells[i].size();j++) {
			index = gCells[i][j];
			if(gFO[index] == MAX_INT) 
				gFO[index] = i;
			else 
				gFN[index] = i;
		}
	}
	/*Flag boundary faces not in gBoundaries for auto deletion*/
	IntVector& gDelete = gBoundaries["delete"];
	for(i = 0;i < gFN.size();i++) {
		if(gFN[i] == MAX_INT) {
			if(!faceInBoundary(i))
				gDelete.push_back(i);
		}
	}
	/*add boundary cells*/
	for(Boundaries::iterator it = gBoundaries.begin();
		it != gBoundaries.end();++it) {
		IntVector& facets = it->second;
		for(j = 0;j < facets.size();j++) {
			i = facets[j];
			//external patch
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
	Int i,j;

	/*allocate*/
	_fC.assign(gFacets.size(),Vector(0));
	_cC.assign(gCells.size(),Vector(0));
	_fN.assign(gFacets.size(),Vector(0));
	_cV.assign(gCells.size(),Scalar(0));
	_reversed.assign(gFacets.size(),false);

	/* face centre*/
	for(i = 0;i < gFacets.size();i++) {
		Facet& f = gFacets[i];
		Vector C(0);
		for(j = 0;j < f.size();j++)
			C += gVertices[f[j]];
		_fC[i] = C / Scalar(f.size());
	}

	/* cell centre */
	for(i = 0;i < gCells.size();i++) {
		Cell& c = gCells[i];
		Vector C(0);
		for(j = 0;j < c.size();j++)
			C += _fC[c[j]];
		_cC[i] = C / Scalar(c.size());
	}
	/* face normal */
	Vector v1,v2,v3,v;
	Scalar magN;
	for(i = 0;i < gFacets.size();i++) {
		Facet& f = gFacets[i];
		Vector N(0),C(0),Ni;
		Scalar Ntot = Scalar(0);
		v1 = _fC[i];
		for(j = 0;j < f.size();j++) {
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
		for(j = 0;j < c.size();j++) {
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

	Int i,j,count;
	IntVector Idf(gFacets.size(),0);
	IntVector Idc(gCells.size(),0);

	/*erase facet reference*/
	for(i = 0;i < fs.size();i++) {
		Int f = fs[i];
		Cell& co = gCells[gFO[f]];
		for(j = 0;j < co.size();j++) {
			if(co[j] == f) {
				co.erase(co.begin() + j); 
				break; 
			}
		}
		Cell& cn = gCells[gFN[f]];
		for(j = 0;j < cn.size();j++) {
			if(cn[j] == f) { 
				cn.erase(cn.begin() + j); 
				break; 
			}
		}
	}
	/*updated facet id*/
	for(i = 0;i < fs.size();i++)
		Idf[fs[i]] = Constants::MAX_INT;
	count = 0;
	for(i = 0;i < gFacets.size();i++) {
		if(Idf[i] != Constants::MAX_INT) 
			Idf[i] = count++;
		else
			gFacets[i].clear();
	}
	/*erase facets*/
	for(i = 0;i < gFacets.size();i++) {
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
	for(i = 0;i < gCells.size();i++) {
		if(gCells[i].size() != 0) 
			Idc[i] = count++;
		else
			Idc[i] = Constants::MAX_INT;
	}
	/*erase cells*/
	for(i = 0;i < gCells.size();i++) {
		if(gCells[i].size() == 0) { 
			gCells.erase(gCells.begin() + i); 
			_cC.erase(_cC.begin() + i);
			_cV.erase(_cV.begin() + i);
			--i;
		} else {
			for(j = 0;j < gCells[i].size();j++) {
				gCells[i][j] = Idf[gCells[i][j]];
			}
		}
	}
	/*facet owner and neighbor*/
	for(i = 0;i < gFacets.size();i++) {
		gFO[i] = Idc[gFO[i]];
		gFN[i] = Idc[gFN[i]];
	}
	/*patches*/
	for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it) {
		IntVector& gB = it->second;
		for(i = 0;i < gB.size();i++) 
			gB[i] = Idf[gB[i]];
	}

	cout << "Total faces: " << gFacets.size() << endl;
}
/*owners*/
IntVector Mesh::owner(const IntVector& f) {
	IntVector r(f.size(),0);
	for(Int i = 0;i < f.size();i++)
		r[i] = gFO[f[i]];
	return r;
}
/*neighbors*/
IntVector Mesh::neighbor(const IntVector& f) {
	IntVector r(f.size(),0);
	for(Int i = 0;i < f.size();i++)
		r[i] = gFN[f[i]];
	return r;
}
/*find nearest vertex*/
int Mesh::findNearest(const Vector& v) {
	Scalar mindist,dist;
	int bi;
	bi = 0;
	mindist = mag(v - gVertices[0]);
	for(Int i = 1;i < gVertices.size();i++) {
		dist = mag(v - gVertices[i]);
		if(dist < mindist) {
			mindist = dist;
			bi = i;
		}
	}
	return bi;
}
