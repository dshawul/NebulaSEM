#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"

namespace Prepare {
	int decompose(Mesh::MeshObject&,Int*,Scalar*,int);
	void decomposeXYZ(Mesh::MeshObject&,Int*,Scalar*,IntVector&);
	void decomposeIndex(Mesh::MeshObject&,Int,IntVector&);
	void decomposeMetis(Mesh::MeshObject&,int,IntVector&);
	void decomposeFields(std::vector<std::string>& fields,std::string,Int);
	int merge(Mesh::MeshObject&,Int*,std::vector<std::string>& fields,std::string,Int);
	int convertVTK(Mesh::MeshObject&,std::vector<std::string>& fields,Int);
	int probe(Mesh::MeshObject&,std::vector<std::string>& fields,Int);
}

#endif
