#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"

namespace Prepare {
	int decompose(Int*,Scalar*,int);
	void decomposeXYZ(Int*,Scalar*,IntVector&);
	void decomposeIndex(Int,IntVector&);
	void decomposeMetis(int,IntVector&);
	void decomposeFields(std::vector<std::string>& fields,std::string,Int);
	int merge(Int*,std::vector<std::string>& fields,std::string,Int);
	int convertVTK(std::vector<std::string>& fields,Int);
	int probe(std::vector<std::string>& fields,Int);
	void refineMesh(const IntVector&,const Int,const Int,const Vector&);
}

#endif
