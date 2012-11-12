#ifndef __DECOMPOSE_H
#define __DECOMPOSE_H

#include "field.h"

namespace Decompose {
	int decomposeXYZ(Mesh::MeshObject&,Int*);
	void decomposeFields(std::vector<std::string>& fields,std::string);
}

#endif
