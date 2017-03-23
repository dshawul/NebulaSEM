#include "tensor.h"

using namespace Constants;

/** Cross product of two vectors */
Vector operator ^ (const Vector& p,const Vector& q) {
    Vector r;
    r[XX] = p[YY] * q[ZZ] - p[ZZ] * q[YY];
    r[YY] = p[ZZ] * q[XX] - p[XX] * q[ZZ];
    r[ZZ] = p[XX] * q[YY] - p[YY] * q[XX];
    return r;
}
/** Outer product of two vectors */
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
/** Outer product of a vector with itself */
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
/** Product of two tensors */
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
/** Product of two symmetric tensors */
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
/** Inner product of a tensor and a vector */
Vector dot(const Tensor& q, const Vector& p) {
    Vector r;
    r[XX] = q[XX] * p[XX] + q[XY] * p[YY] + q[XZ] * p[ZZ];
    r[YY] = q[YX] * p[XX] + q[YY] * p[YY] + q[YZ] * p[ZZ];
    r[ZZ] = q[ZX] * p[XX] + q[ZY] * p[YY] + q[ZZ] * p[ZZ];
    return r;
}
/** Inner product of a symmetric tensor and a vector */
Vector dot(const STensor& q, const Vector& p) {
    Vector r;
    r[XX] = q[XX] * p[XX] + q[XY] * p[YY] + q[XZ] * p[ZZ];
    r[YY] = q[XY] * p[XX] + q[YY] * p[YY] + q[YZ] * p[ZZ];
    r[ZZ] = q[XZ] * p[XX] + q[YZ] * p[YY] + q[ZZ] * p[ZZ];
    return r;
}
/** Symmetric part of a tensor */
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
/** Skew-symmetric part of a tensor */
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
/** Transpose of a tensor */
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
/** Rotate a vector by an angle */
Vector rotate(const Vector& v,const Vector& N,const Scalar& theta) {
    Vector r;
    Scalar sum = v & N;
    Scalar cost = cos(theta), sint = sin(theta);
    r[XX] = N[XX] * sum * (1 - cost) + v[XX] * cost + (-N[ZZ] * v[YY] + N[YY] * v[ZZ]) * sint;
    r[YY] = N[YY] * sum * (1 - cost) + v[YY] * cost + (+N[ZZ] * v[XX] - N[XX] * v[ZZ]) * sint;
    r[ZZ] = N[ZZ] * sum * (1 - cost) + v[ZZ] * cost + (-N[YY] * v[XX] + N[XX] * v[YY]) * sint;
    return r;
}
/** Determinant of a tensor */
Scalar det(const Tensor& p) {
    Scalar r;
    r =  p[XX] * (p[YY] * p[ZZ] - p[YZ] * p[ZY]) +
         p[XY] * (p[YZ] * p[ZX] - p[YX] * p[ZZ]) +
         p[XZ] * (p[YX] * p[ZY] - p[YY] * p[ZX]);
    return r;
}
/** Inverse of a tensor */
Tensor inv(const Tensor& p) {
    Tensor r;
    r[XX] = p[YY] * p[ZZ] - p[YZ] * p[ZY];
    r[YY] = p[XX] * p[ZZ] - p[XZ] * p[ZX];
    r[ZZ] = p[XX] * p[YY] - p[XY] * p[YX];
    r[XY] = p[XZ] * p[ZY] - p[XY] * p[ZZ];
    r[XZ] = p[XY] * p[YZ] - p[XZ] * p[YY];
    r[YX] = p[YZ] * p[ZX] - p[YX] * p[ZZ];
    r[YZ] = p[XZ] * p[YX] - p[XX] * p[YZ];
    r[ZX] = p[YX] * p[ZY] - p[YY] * p[ZX];
    r[ZY] = p[XY] * p[ZX] - p[XX] * p[ZY];
    Scalar d = p[XX] * r[XX] + p[XY] * r[YX] + p[XZ] * r[ZX];
    if(d == 0) r = Tensor(0);
    else r /= d;
    return r;
}
