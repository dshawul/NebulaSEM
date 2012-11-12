#include "tensor.h"

using namespace Constants;

Vector operator ^ (const Vector& p,const Vector& q) {
	Vector r;
	r[XX] = p[YY] * q[ZZ] - p[ZZ] * q[YY];
	r[YY] = p[ZZ] * q[XX] - p[XX] * q[ZZ];
	r[ZZ] = p[XX] * q[YY] - p[YY] * q[XX];
	return r;
}
Tensor mul(const Vector& p,const Vector& q) {
	Tensor r;
	r[XX] = p[XX] * q[XX];
	r[YY] = p[YY] * q[YY];
	r[ZZ] = p[ZZ] * q[ZZ];

	r[XY] = p[XX] * q[YY];
	r[YZ] = p[YY] * q[ZZ];
	r[XZ] = p[XX] * q[ZZ];

	r[YX] = p[YY] * q[XX];
	r[ZY] = p[ZZ] * q[YY];
	r[ZX] = p[ZZ] * q[XX];
	return r;
}
STensor mul(const Vector& p) {
	STensor r;
	r[XX] = p[XX] * p[XX];
	r[YY] = p[YY] * p[YY];
	r[ZZ] = p[ZZ] * p[ZZ];

	r[XY] = p[XX] * p[YY];
	r[YZ] = p[YY] * p[ZZ];
	r[XZ] = p[XX] * p[ZZ];
	return r;
}
Tensor mul(const Tensor& p,const Tensor& q) {
	Tensor r;
	r[XX] = p[XX] * q[XX] + p[XY] * q[YX] + p[XZ] * q[ZX];
	r[XY] = p[XX] * q[XY] + p[XY] * q[YY] + p[XZ] * q[ZY];
	r[XZ] = p[XX] * q[XZ] + p[XY] * q[YZ] + p[XZ] * q[ZZ];

	r[YX] = p[YX] * q[XX] + p[YY] * q[YX] + p[YZ] * q[ZX];
	r[YY] = p[YX] * q[XY] + p[YY] * q[YY] + p[YZ] * q[ZY];
	r[YZ] = p[YX] * q[XZ] + p[YY] * q[YZ] + p[YZ] * q[ZZ];

	r[ZX] = p[ZX] * q[XX] + p[ZY] * q[YX] + p[ZZ] * q[ZX];
	r[ZY] = p[ZX] * q[XY] + p[ZY] * q[YY] + p[ZZ] * q[ZY];
	r[ZZ] = p[ZX] * q[XZ] + p[ZY] * q[YZ] + p[ZZ] * q[ZZ];

	return r;
}
STensor mul(const STensor& p,const STensor& q) {
	STensor r;
	r[XX] = p[XX] * q[XX] + p[XY] * q[XY] + p[XZ] * q[XZ];
	r[XY] = p[XX] * q[XY] + p[XY] * q[YY] + p[XZ] * q[YZ];
	r[XZ] = p[XX] * q[XZ] + p[XY] * q[YZ] + p[XZ] * q[ZZ];

	r[YY] = p[XY] * q[XY] + p[YY] * q[YY] + p[YZ] * q[YZ];
	r[YZ] = p[XY] * q[XZ] + p[YY] * q[YZ] + p[YZ] * q[ZZ];
	r[ZZ] = p[XZ] * q[XZ] + p[YZ] * q[YZ] + p[ZZ] * q[ZZ];

	return r;
}
Vector dot(const Vector& p,const Tensor& q) {
	Vector r;
	r[XX] = q[XX] * p[XX] + q[XY] * p[YY] + q[XZ] * p[ZZ];
	r[YY] = q[YX] * p[XX] + q[YY] * p[YY] + q[YZ] * p[ZZ];
	r[ZZ] = q[ZX] * p[XX] + q[ZY] * p[YY] + q[ZZ] * p[ZZ];
	return r;
}
Vector dot(const Vector& p,const STensor& q) {
	Vector r;
	r[XX] = q[XX] * p[XX] + q[XY] * p[YY] + q[XZ] * p[ZZ];
	r[YY] = q[XY] * p[XX] + q[YY] * p[YY] + q[YZ] * p[ZZ];
	r[ZZ] = q[XZ] * p[XX] + q[YZ] * p[YY] + q[ZZ] * p[ZZ];
	return r;
}
STensor sym(const Tensor& p) {
	STensor r;
	r[XX] = p[XX];
	r[YY] = p[YY];
	r[ZZ] = p[ZZ];

	r[XY] = (p[XY] + p[YX]) / 2;
	r[YZ] = (p[YZ] + p[ZY]) / 2;
	r[XZ] = (p[XZ] + p[ZX]) / 2;
	return r;
}
Tensor skw(const Tensor& p) {
	Tensor r;
	r[XX] = 0;
	r[YY] = 0;
	r[ZZ] = 0;

	r[XY] = (p[XY] - p[YX]) / 2;
	r[YZ] = (p[YZ] - p[ZY]) / 2;
	r[XZ] = (p[XZ] - p[ZX]) / 2;

	r[YX] = (p[YX] - p[XY]) / 2;
	r[ZY] = (p[ZY] - p[YZ]) / 2;
	r[ZX] = (p[ZX] - p[XZ]) / 2;
	return r;
}
Tensor trn(const Tensor& p) {
	Tensor r = p;
	r[XX] = p[XX];
	r[YY] = p[YY];
	r[ZZ] = p[ZZ];

	r[XY] = p[YX];
	r[YZ] = p[ZY];
	r[XZ] = p[ZX];

	r[YX] = p[XY];
	r[ZY] = p[YZ];
	r[ZX] = p[XZ];
	return r;
}
