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

	/*add boundary cells*/
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
	/* 1. Add boundary cells
	 * 2. Faces not in gBoundaries are flagged for auto deletion
	 */
	IntVector& gDelete = gBoundaries["delete"];
	for(i = 0;i < gFN.size();i++) {
		if(gFN[i] == MAX_INT) {
			/*add to delete list*/
			if(!faceInBoundary(i))
				gDelete.push_back(i);
			/*add it*/
			Cell c;
			c.push_back(i);
			gCells.push_back(c);
			gFN[i] = gCells.size() - 1;
		}
	}
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
