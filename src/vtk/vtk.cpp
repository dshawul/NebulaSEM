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
    Int i,total;
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
    
    if(DG::NPMAT <= 1) {
        /*vertices*/
        {
            of << "POINTS " << gVertices.size() << " double" << endl;
            of.precision(12);
            forEach(gVertices,i)
                of << gVertices[i] << endl;
            of.precision(6);
        }
        /*cells*/
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
            for(i = 0;i < gBCS;i++)
                of << "12" << endl;
        }
        /*Fields*/
        total = ScalarCellField::count_writable() +
                VectorCellField::count_writable() +
                STensorCellField::count_writable() +
                TensorCellField::count_writable();
            
        /*cells*/
        if(write_cell_value) {
            of << "CELL_DATA " << gBCS << endl;
            of << "FIELD attributes "<< total + 1 << endl;
            forEachCellField(writeVtkCellAll(of));
            of << "cellID  1 " << gBCS << " int" << endl;
            for(Int i = 0;i < gBCS;i++) of << i << endl;
        }
        /*vertices*/
        {
            of << "POINT_DATA " << gVertices.size() << endl;
            of << "FIELD attributes "<< total << endl;
            forEachCellField(writeVtkVertexAll(of));
        }
    } else {
        using namespace DG;
        /*vertices*/
        {
            of << "POINTS " << gBCSfield << " double" << endl;
            of.precision(12);
            for(i = 0;i < gBCSfield;i++)
                of << cC[i] << endl;
            of.precision(6);
        }
        /*points,lines, polygons and hexahedras*/
        {
            //points
            if(NPX == 1 && NPY == 1 && NPZ == 1) {
                Int ncells = gBCS;
                of << "CELLS " << ncells << " " << ncells * 2 << endl;
                Int r = 0, s = 0, t = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                {
                    of << "1 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "1" << endl;
            //lines
            } else if(NPX == 1 && NPY == 1) {
                Int ncells = gBCS * (NPZ - 1);
                of << "CELLS " << ncells << " " << ncells * 3 << endl;
                Int r = 0, s = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int t = 0;t < NPZ - 1;t++)
                {
                    of << "2 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,  r,  s,  t+1) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "3" << endl;
            } else if(NPX == 1 && NPZ == 1) {
                Int ncells = gBCS * (NPY - 1);
                of << "CELLS " << ncells << " " << ncells * 3 << endl;
                Int r = 0, t = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int s = 0;s < NPY - 1;s++)
                {
                    of << "2 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,  r,s+1,  t) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "3" << endl;
            } else  if(NPY == 1 && NPZ == 1) {
                Int ncells = gBCS * (NPX - 1);
                of << "CELLS " << ncells << " " << ncells * 3 << endl;
                Int s = 0, t = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int r = 0;r < NPX - 1;r++)
                {
                    of << "2 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,r+1,  s,  t) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "3" << endl;
            //polygons
            } else if(NPX == 1) {
                Int ncells = gBCS * (NPY - 1) * (NPZ - 1);
                of << "CELLS " << ncells << " " << ncells * 5 << endl;
                Int r = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int s = 0;s < NPY - 1;s++)
                for(Int t = 0;t < NPZ - 1;t++)
                {
                    of << "4 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,  r,s+1,  t) << " ";
                    of << INDEX4(ci,  r,s+1,  t+1) << " ";
                    of << INDEX4(ci,  r,  s,  t+1) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "9" << endl;
            } else if(NPY == 1) {
                Int ncells = gBCS * (NPX - 1) * (NPZ - 1);
                of << "CELLS " << ncells << " " << ncells * 5 << endl;
                Int s = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int r = 0;r < NPX - 1;r++)
                for(Int t = 0;t < NPZ - 1;t++)
                {
                    of << "4 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,r+1,  s,  t) << " ";
                    of << INDEX4(ci,r+1,  s,  t+1) << " ";
                    of << INDEX4(ci,  r,  s,  t+1) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "9" << endl;
            } else if(NPZ == 1) {
                Int ncells = gBCS * (NPX - 1) * (NPY - 1);
                of << "CELLS " << ncells << " " << ncells * 5 << endl;
                Int t = 0;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int r = 0;r < NPX - 1;r++)
                for(Int s = 0;s < NPY - 1;s++)
                {
                    of << "4 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,r+1,  s,  t) << " ";
                    of << INDEX4(ci,r+1,s+1,  t) << " ";
                    of << INDEX4(ci,  r,s+1,  t) << " ";
                    of << endl; 
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "9" << endl;
            //cells
            } else {
                Int ncells = gBCS * (NPX - 1) * (NPY - 1) * (NPZ - 1);
                of << "CELLS " << ncells << " " << ncells * 9 << endl;
                for(Int ci = 0;ci < gBCS;ci++)
                for(Int r = 0;r < NPX - 1;r++)
                for(Int s = 0;s < NPY - 1;s++)
                for(Int t = 0;t < NPZ - 1;t++)
                {
                    of << "8 ";
                    of << INDEX4(ci,  r,  s,  t) << " ";
                    of << INDEX4(ci,r+1,  s,  t) << " ";
                    of << INDEX4(ci,r+1,s+1,  t) << " ";
                    of << INDEX4(ci,  r,s+1,  t) << " ";
                    of << INDEX4(ci,  r,  s,t+1) << " ";
                    of << INDEX4(ci,r+1,  s,t+1) << " ";
                    of << INDEX4(ci,r+1,s+1,t+1) << " ";
                    of << INDEX4(ci,  r,s+1,t+1) << " ";
                    of << endl;
                }
                of << "CELL_TYPES " << ncells << endl;
                for(i = 0;i < ncells;i++)
                    of << "12" << endl;
            }
        }
        /*Fields*/
        total = ScalarCellField::count_writable() +
                VectorCellField::count_writable() +
                STensorCellField::count_writable() +
                TensorCellField::count_writable();
        /*vertices*/
        {
            of << "POINT_DATA " << gBCSfield << endl;
            of << "FIELD attributes "<< total << endl;
            forEachCellField(writeVtkCellAll(of));
        }
    }
}
