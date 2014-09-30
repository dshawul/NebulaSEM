#ifndef __DG_H
#define __DG_H

#include "tensor.h"

namespace DG {
	void legendre(int p, Scalar x,Scalar& L0,Scalar& L0_1,Scalar& L0_2);
	void LGL(int N, Scalar* xgl, Scalar* wgl);
}

#endif