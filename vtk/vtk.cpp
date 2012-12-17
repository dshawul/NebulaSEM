#include "vtk.h"

using namespace std;
using namespace Mesh;

//#define POLYHEDRAL

#ifdef POLYHEDRAL
static Int cell_count(Cell& c) {
	Facet* f;
	Int i,nFacets = c.size(),nVertices = 0,nTotal;
	for(i = 0;i < nFacets;i++) {
		f = &gFacets[c[i]];
		nVertices += f->size();
	}
	nTotal = nFacets + nVertices + 2;
	return nTotal;
}

static void cell_vtk(std::ofstream& of, Cell& c) {
	Facet* f;
	Int i,j,nFacets = c.size(),nVertices = 0,nTotal;
	for(i = 0;i < nFacets;i++) {
		f = &gFacets[c[i]];
		nVertices += f->size();
	}
	nTotal = nFacets + nVertices + 2;

	/*write*/
	of << nTotal - 1 << " " << nFacets << " ";
	for(i = 0;i < nFacets;i++) {
		f = &gFacets[c[i]];
		of << f->size() << " ";
		for(j = 0;j < f->size();j++) {
			of << (*f)[j] << " ";
		}
	}
	of << endl;
}
#endif

void Util::write_vtk(Int step) {
	Int total;
	stringstream path;
	path << gMeshName << step << ".vtk";
	ofstream of(path.str().c_str());
#ifdef POLYHEDRAL
	of << "# vtk DataFile Version 2.0" << endl;
#else
	of << "# vtk DataFile Version 1.0" << endl;
#endif
	of << Mesh::gMeshName << endl;
	of << "ASCII" << endl;
	of << "DATASET UNSTRUCTURED_GRID" << endl;
    /*Geometry*/
	Int i,j;
	of << "POINTS " << gVertices.size() << " float" << endl;
	for(i = 0;i < gVertices.size();i++) {
		of << gVertices[i] << endl;
    }
#ifdef POLYHEDRAL
	/*polyhedral cells*/
	total = 0;
	for(i = 0;i < gBCellsStart;i++)
		total += cell_count(gCells[i]);

	of << "CELLS " << gBCellsStart << " " << total << endl;
	for(i = 0;i < gBCellsStart;i++)
		cell_vtk(of,gCells[i]);
    
	of << "CELL_TYPES " << gBCellsStart << endl;
	for(i = 0;i < gBCellsStart;i++)
		of << 42 << endl;
#else
    /*hexahedral cells*/
    of << "CELLS " << gBCellsStart << " " << gBCellsStart * 9 << endl;
	for(i = 0;i < gBCellsStart;i++) {
		Cell& c = gCells[i];
		Facet f1 = gFacets[c[0]];
		Facet f2 = gFacets[c[1]];
		of << f1.size() + f2.size() << " ";
		for(j = 0;j < f1.size();j++) of << f1[j] << " ";
		for(j = 0;j < f2.size();j++) of << f2[j] << " ";
		of << endl;
    }
	of << "CELL_TYPES " << gBCellsStart << endl;
	for(i = 0;i < gBCellsStart;i++) {
		of << "12" << endl;
    }
#endif
	/*Fields*/
	total = ScalarCellField::count_writable() +
		    VectorCellField::count_writable() +
			STensorCellField::count_writable() +
			TensorCellField::count_writable();

	of << "CELL_DATA " << gBCellsStart << endl;
	of << "FIELD attributes "<< total + 1 << endl;
	ScalarCellField::write_vtk(of,false);
	VectorCellField::write_vtk(of,false);
	STensorCellField::write_vtk(of,false);
	TensorCellField::write_vtk(of,false);
	of << "cellID  1 " << Mesh::gBCellsStart << " int" << endl;
	for(Int i = 0;i < Mesh::gBCellsStart;i++) of << i << endl;

	of << "POINT_DATA " << gVertices.size() << endl;
	of << "FIELD attributes "<< total << endl;
	ScalarCellField::write_vtk(of,true);
	VectorCellField::write_vtk(of,true);
	STensorCellField::write_vtk(of,true);
	TensorCellField::write_vtk(of,true);

	/*end*/
	of.close();
}
