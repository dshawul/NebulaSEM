#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"

namespace Prepare {
	int decompose(Mesh::MeshObject&,Int*,Scalar*);
	void decomposeXYZ(Mesh::MeshObject&,Int*,Scalar*,IntVector&);
	void decomposeFields(std::vector<std::string>& fields,std::string,Int);
	int merge(Mesh::MeshObject&,Int*,std::vector<std::string>& fields,std::string,Int);
	int convertVTK(Mesh::MeshObject&,std::vector<std::string>& fields,Int);
	int probe(Mesh::MeshObject&,std::vector<std::string>& fields,Int);
}

#endif
