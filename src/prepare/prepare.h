#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"
	
namespace Prepare {
	void decomposeXYZ(Int*,Scalar*,IntVector&);
	void decomposeIndex(Int,IntVector&);
	void decomposeMetis(int,IntVector&);
	int decompose(std::vector<std::string>&,Int*,Scalar*,int,Int);
	int merge(Int*,std::vector<std::string>&,Int);
	int convertVTK(std::vector<std::string>&,Int);
	int probe(std::vector<std::string>&,Int);
}

#endif
