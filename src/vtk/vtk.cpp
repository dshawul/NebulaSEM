#include "vtk.h"

using namespace std;
using namespace Mesh;

bool Vtk::write_polyhedral = false;
bool Vtk::write_cell_value = true;

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
	Int i,j,nFacets = c.size(),nTotal = cell_count(c);
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

void Vtk::write_vtk(Int step) {
	Int total;
	stringstream path;
	path << gMeshName << step << ".vtk";
	ofstream of(path.str().c_str());
	if(write_polyhedral)
		of << "# vtk DataFile Version 2.0" << endl;
	else
		of << "# vtk DataFile Version 1.0" << endl;
	of << Mesh::gMeshName << endl;
	of << "ASCII" << endl;
	of << "DATASET UNSTRUCTURED_GRID" << endl;
    /*Geometry*/
	Int i;
	of << "POINTS " << gVertices.size() << " float" << endl;
	of.precision(12);
	forEach(gVertices,i)
		of << gVertices[i] << endl;
    of.precision(6);
	if(write_polyhedral) {
		/*polyhedral cells*/
		total = 0;
		for(i = 0;i < gBCS;i++)
			total += cell_count(gCells[i]);

		of << "CELLS " << gBCS << " " << total << endl;
		for(i = 0;i < gBCS;i++)
			cell_vtk(of,gCells[i]);

		of << "CELL_TYPES " << gBCS << endl;
		for(i = 0;i < gBCS;i++)
			of << "42" << endl;
	} else {
		/*hexahedral cells*/
		of << "CELLS " << gBCS << " " << gBCS * 9 << endl;
		for(i = 0;i < gBCS;i++) {
			Cell& c = gCells[i];
			Facet f1 = gFacets[c[0]];
			Facet f2 = gFacets[c[1]];
			of << f1.size() + f2.size() << " ";
			forEach(f1,j) 
				of << f1[j] << " ";
			forEach(f2,j)
				of << f2[j] << " ";
			of << endl;
		}
		of << "CELL_TYPES " << gBCS << endl;
		for(i = 0;i < gBCS;i++) {
			of << "12" << endl;
		}
	}
	/*Fields*/
	total = ScalarCellField::count_writable() +
		    VectorCellField::count_writable() +
			STensorCellField::count_writable() +
			TensorCellField::count_writable();
	if(write_cell_value) {
		of << "CELL_DATA " << gBCS << endl;
		of << "FIELD attributes "<< total + 1 << endl;
		forEachCellField(writeVtkCellAll(of));
		of << "cellID  1 " << gBCS << " int" << endl;
		for(Int i = 0;i < gBCS;i++) of << i << endl;
	}
	of << "POINT_DATA " << gVertices.size() << endl;
	of << "FIELD attributes "<< total << endl;
	forEachCellField(writeVtkVertexAll(of));
}
