#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"

struct RefineParams {
	Int shape;
	Vector dir;
	std::string field;
	Scalar field_max;
	Scalar field_min;
	Int limit;
	RefineParams() {
		shape = 0;
		dir = Scalar(0);
		field = "U";
		field_max = 0.9;
		field_min = 0.1;
		limit = 100000;
	}
};
	
namespace Prepare {
	int decompose(std::vector<std::string>&,Int*,Scalar*,int,Int);
	void decomposeXYZ(Int*,Scalar*,IntVector&);
	void decomposeIndex(Int,IntVector&);
	void decomposeMetis(int,IntVector&);
	int merge(Int*,std::vector<std::string>& fields,Int);
	int convertVTK(std::vector<std::string>& fields,Int);
	int probe(std::vector<std::string>&,Int);
	void refineMesh(std::vector<std::string>&,const RefineParams&,Int);
	void refineMesh(IntVector&,const RefineParams&,IntVector&);
}

#endif
