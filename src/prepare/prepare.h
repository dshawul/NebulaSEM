#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"
	
namespace Prepare {
	int decompose(std::vector<std::string>&,Int*,Scalar*,int,Int);
	void decomposeXYZ(Int*,Scalar*,IntVector&);
	void decomposeIndex(Int,IntVector&);
	void decomposeMetis(int,IntVector&);
	int merge(Int*,std::vector<std::string>& fields,Int);
	int convertVTK(std::vector<std::string>& fields,Int);
	int probe(std::vector<std::string>&,Int);
	void refineMesh(std::vector<std::string>&,const RefineParams&,Int);
}

#endif
