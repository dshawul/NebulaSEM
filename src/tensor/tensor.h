#ifndef __TENSOR_H
#define __TENSOR_H

#include <iostream>
#include <cmath>

#include "types.h"

/**************************
 * scalars
 *************************/
#if defined USE_DOUBLE
#   define Scalar  double
#else
#   define Scalar  float
#endif

/** Adds general arthimetic operators defined via compound assignment*/
#define AddOperator(T,$)                                            \
    friend FORCEINLINE T operator $ (const T& p,const T& q) {       \
        T r = p;                                                    \
        r $##= q;                                                   \
        return r;                                                   \
    }
/** \name General scalar operators */
//@{
#define AddRightScalarOperator(T,$)                                 \
    friend FORCEINLINE T operator $ (const T& p,const Scalar& q) {  \
        T r = p;                                                    \
        r $##= q;                                                   \
        return r;                                                   \
    }
#define AddLeftScalarOperator1(T,$)                                 \
    friend FORCEINLINE T operator $ (const Scalar& p,const T& q) {  \
        T r = q;                                                    \
        r $##= p;                                                   \
        return r;                                                   \
    }
#define AddLeftScalarOperator2(T,$)                                 \
    friend FORCEINLINE T operator $ (const Scalar& p,const T& q) {  \
        T r(p);                                                     \
        r $##= q;                                                   \
        return r;                                                   \
    }
//@}
    
#define AddOperators(T)\
    AddOperator(T,+)\
    AddOperator(T,-)\
    AddOperator(T,*)\
    AddOperator(T,/)
    
#define AddScalarOperators(T)\
    AddRightScalarOperator(T,+)\
    AddRightScalarOperator(T,-)\
    AddRightScalarOperator(T,*)\
    AddRightScalarOperator(T,/)\
    AddLeftScalarOperator1(T,+)\
    AddLeftScalarOperator1(T,*)\
    AddLeftScalarOperator2(T,-)\
    AddLeftScalarOperator2(T,/)

/** \name Other operations on scalars*/
//@{
FORCEINLINE Scalar* addr(Int& p) { 
    return (Scalar*) &p; 
}
FORCEINLINE Scalar* addr(Scalar& p) { 
    return &p; 
}
FORCEINLINE Scalar mag(const Scalar& p) { 
    return fabs(p); 
}
FORCEINLINE Scalar sdiv(const Scalar& p,const Scalar& q) { 
    return p ? (p / q) : 0; 
}
FORCEINLINE Scalar max(const Scalar& p,const Scalar& q) { 
    return (p >= q) ? p : q; 
}
FORCEINLINE Scalar min(const Scalar& p,const Scalar& q) { 
    return (p <= q) ? p : q; 
}
FORCEINLINE Scalar dot(const Scalar& p, const Scalar& q) {
    return p * q;
}
//@}

/*********************************
 * Tensors
 *********************************/

/**
Loop unroller for the Tensor class
*/
template <int N>
struct Unroll {
    /*macro*/
#define Op(name,$)                                                              \
    static FORCEINLINE void name(Scalar* p,const Scalar* q) {                   \
        *p $ *q;                                                                \
        Unroll<N - 1>::name(p + 1,q + 1);                                       \
    }     
#define SOp(name,$)                                                             \
    static FORCEINLINE void name(Scalar* p,const Scalar q) {                    \
        *p $ q;                                                                 \
        Unroll<N - 1>::name(p + 1,q);                                           \
    }
#define Fp(name,$)                                                              \
    static FORCEINLINE void name(Scalar* r,const Scalar* p,const Scalar* q) {   \
        *r = ::$(*p,*q);                                                        \
        Unroll<N - 1>::name(r + 1,p + 1,q + 1);                                 \
    }
#define Fp1(name,$)                                                             \
    static FORCEINLINE void name(Scalar* r,const Scalar* p,const Scalar q) {    \
        *r = ::$(*p,q);                                                         \
        Unroll<N - 1>::name(r + 1,p + 1,q);                                     \
    }
#define Fp2(name,$)                                                             \
    static FORCEINLINE void name(Scalar* r,const Scalar* p) {                   \
        *r = ::$(*p);                                                           \
        Unroll<N - 1>::name(r + 1,p + 1);                                       \
    }
    /*special*/
    static FORCEINLINE Scalar dot(const Scalar* p,const Scalar* q) {     
        return (*p) * (*q) + Unroll<N - 1>::dot(p + 1,q + 1);                                    
    }
    /*define ops*/
    Op(equ,=)
    Op(neg,=-)
    Op(inc,+=)
    Op(dec,-=)
    Op(mul,*=)
    Op(div,/=)
    SOp(equ,=)
    SOp(neg,=-)
    SOp(inc,+=)
    SOp(dec,-=)
    SOp(mul,*=)
    SOp(div,/=)
    Fp(sdiv,sdiv)
    /*from math.h*/
    Fp2(acos,acos)
    Fp2(asin,asin)
    Fp2(atan,atan)
    Fp(atan2,atan2)
    Fp2(ceil,ceil)
    Fp2(cos,cos)
    Fp2(cosh,cosh)
    Fp2(exp,exp)
    Fp2(fabs,fabs)
    Fp2(floor,floor)
    Fp2(log,log)
    Fp2(log10,log10)
    Fp1(pow,pow)
    Fp2(sin,sin)
    Fp2(sinh,sinh)
    Fp2(sqrt,sqrt)
    Fp2(tan,tan)
    Fp2(tanh,tanh)
    Fp(min,min)
    Fp(max,max)
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
};

/**
Terminator class for the loop unroller
*/
template <>
struct Unroll<0> {
    /*macro*/
#define Op(name)                                                              \
    static FORCEINLINE void name(Scalar* p,const Scalar* q) {}   
#define SOp(name)                                                             \
    static FORCEINLINE void name(Scalar* p,const Scalar q) {}   
#define Fp(name)                                                              \
    static FORCEINLINE void name(Scalar* r,const Scalar* p,const Scalar* q) {}   
#define Fp1(name)                                                             \
    static FORCEINLINE void name(Scalar* r,const Scalar* p,const Scalar q) {}   
#define Fp2(name)                                                             \
    static FORCEINLINE void name(Scalar* r,const Scalar* p) {}   
    /*special*/
    static FORCEINLINE Scalar dot(const Scalar* p,const Scalar* q) {return 0;}
    /*define ops*/
    Op(equ)
    Op(neg)
    Op(inc)
    Op(dec)
    Op(mul)
    Op(div)
    SOp(equ)
    SOp(neg)
    SOp(inc)
    SOp(dec)
    SOp(mul)
    SOp(div)
    Fp(sdiv)
    /*from math.h*/
    Fp2(acos)
    Fp2(asin)
    Fp2(atan)
    Fp(atan2)
    Fp2(ceil)
    Fp2(cos)
    Fp2(cosh)
    Fp2(exp)
    Fp2(fabs)
    Fp2(floor)
    Fp2(log)
    Fp2(log10)
    Fp1(pow)
    Fp2(sin)
    Fp2(sinh)
    Fp2(sqrt)
    Fp2(tan)
    Fp2(tanh)
    Fp(min)
    Fp(max)
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
};


/**
Template tensor class for defining different size
tensors Scalar=1 Vector=3, Symmetric tensor = 6 and
asymmetric tensor = 9
*/
template <Int SIZE>
class TTensor {
public:
    Scalar P[SIZE];
public:
    /*c'tors*/
    TTensor() {
    }
    TTensor(const TTensor& p) {
        *this = p;
    }
    explicit TTensor(const Scalar& p) {
        Unroll<SIZE>::equ(P,p);
    }
    TTensor(const Scalar& xx,const Scalar& yy,const Scalar& zz) {
        P[0] = xx;
        P[1] = yy;
        P[2] = zz;
        Unroll<SIZE - 3>::equ(&P[3],Scalar(0));
    }
    /*accessors*/
    FORCEINLINE Scalar& operator [] (Int i) {
        return P[i];
    }
    FORCEINLINE const Scalar& operator [] (Int i) const {
        return P[i];
    }
    /*unary ops*/
    FORCEINLINE TTensor operator - () {
        TTensor r;
        Unroll<SIZE>::neg(r.P,P);
        return r;
    }
    FORCEINLINE friend Scalar operator & (const TTensor& p,const TTensor& q) {
        Scalar r = Unroll<SIZE>::dot(p.P,q.P);
        if(SIZE == 6) r += Unroll<3>::dot(&p.P[3],&q.P[3]);
        return r;
    }
    FORCEINLINE friend Scalar* addr(TTensor& p) {
        return p.P;
    } 
    /*unrolled operations*/
#define Op(name,$)                                                  \
    FORCEINLINE TTensor& operator $(const TTensor& q) {             \
        Unroll<SIZE>::name(P,q.P);                                  \
        return *this;                                               \
    }
#define SOp(name,$)                                                 \
    FORCEINLINE TTensor& operator $(const Scalar& q) {              \
        Unroll<SIZE>::name(P,q);                                    \
        return *this;                                               \
    }
#define Fp(name)                                                    \
    FORCEINLINE friend TTensor name(const TTensor& p,const TTensor& s) { \
        TTensor r;                                                  \
        Unroll<SIZE>::name(r.P,p.P,s.P);                            \
        return r;                                                   \
    }
#define Fp1(name)                                                   \
    FORCEINLINE friend TTensor name(const TTensor& p,const Scalar& s) { \
        TTensor r;                                                  \
        Unroll<SIZE>::name(r.P,p.P,s);                              \
        return r;                                                   \
    }
#define Fp2(name)                                                   \
    FORCEINLINE friend TTensor name(const TTensor& p) {             \
        TTensor r;                                                  \
        Unroll<SIZE>::name(r.P,p.P);                                \
        return r;                                                   \
    }
    /*define ops*/
    Op(equ,=)
    Op(inc,+=)
    Op(dec,-=)
    Op(mul,*=)
    Op(div,/=)
    SOp(equ,=)
    SOp(inc,+=)
    SOp(dec,-=)
    SOp(mul,*=)
    SOp(div,/=)
    Fp(sdiv)
    /*from math.h*/
    Fp2(acos)
    Fp2(asin)
    Fp2(atan)
    Fp(atan2)
    Fp2(ceil)
    Fp2(cos)
    Fp2(cosh)
    Fp2(exp)
    Fp2(fabs)
    Fp2(floor)
    Fp2(log)
    Fp2(log10)
    Fp1(pow)
    Fp2(sin)
    Fp2(sinh)
    Fp2(sqrt)
    Fp2(tan)
    Fp2(tanh)
    Fp(min)
    Fp(max)
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
    AddOperators(TTensor)
    AddScalarOperators(TTensor)
    /*others*/
    FORCEINLINE friend Scalar dot(const TTensor& p,const TTensor& q) {
        return p & q;
    }
    FORCEINLINE friend TTensor mul(const TTensor& p,const Scalar& q) {
        return p * q;
    }
    FORCEINLINE friend Scalar magSq(const TTensor& p) {
        return p & p;
    }
    FORCEINLINE friend Scalar mag(const TTensor& p) {
        return sqrt(magSq(p));
    }
    FORCEINLINE friend TTensor unit(const TTensor& p) {
        TTensor r = p;
        Scalar mg = mag(r);
        r /= mg;
        return r;
    }
    FORCEINLINE friend Scalar tr(const TTensor& p) {
        return p[0] + p[1] + p[2];
    }
    FORCEINLINE friend TTensor dev(const TTensor& p,const Scalar& factor) {
        TTensor r = p;
        Scalar t = tr(p) * factor / 3;
        r[0] -= t;
        r[1] -= t;
        r[2] -= t;
        return r;
    }
    FORCEINLINE friend TTensor hyd(const TTensor& p,const Scalar& factor) {
        TTensor r(1,1,1);
        Scalar t = tr(p) * factor / 3;
        r[0] = t;
        r[1] = t;
        r[2] = t;
        return r;
    }
    /*IO*/
    template<typename Ts>
    friend Ts& operator << (Ts& os, const TTensor<SIZE>& p) {
        for(Int i = 0;i < SIZE;i++) 
            os << p[i] << " ";
        return os;
    }
    template<typename Ts>
    friend Ts& operator >> (Ts& is, TTensor<SIZE>& p) {
        for(Int i = 0;i < SIZE;i++) 
            is >> p[i];
        return is;
    }
};

/** \name Define basic tensors of order 3,6, and 9*/
//@{
typedef TTensor<3> Vector;
typedef TTensor<6> STensor;
typedef TTensor<9> Tensor;
//@}

/** OpenMP sum reductions */
#pragma omp declare	\
reduction(+ : Vector : omp_out += omp_in )	\
initializer( omp_priv = omp_orig )

#pragma omp declare	\
reduction(+ : STensor : omp_out += omp_in )	\
initializer( omp_priv = omp_orig )

#pragma omp declare	\
reduction(+ : Tensor : omp_out += omp_in )	\
initializer( omp_priv = omp_orig )

/*Tensor operations*/
#pragma acc routine seq
Vector operator ^ (const Vector& p,const Vector& q);
#pragma acc routine seq
STensor mul(const Vector& p);
#pragma acc routine seq
Tensor mul(const Vector& p,const Vector& q);
#pragma acc routine seq
Tensor mul(const Tensor& p,const Tensor& q);
#pragma acc routine seq
STensor mul(const STensor& p,const STensor& q);
#pragma acc routine seq
Vector dot(const Tensor&,const Vector&);
#pragma acc routine seq
Vector dot(const STensor&,const Vector&);
#pragma acc routine seq
STensor sym(const Tensor& p);
#pragma acc routine seq
Tensor skw(const Tensor& p);
#pragma acc routine seq
Tensor trn(const Tensor& p);
#pragma acc routine seq
Scalar det(const Tensor& p);
#pragma acc routine seq
Tensor inv(const Tensor& p);
#pragma acc routine seq
Vector rotate(const Vector& v,const Vector& N,const Scalar& theta); 

/**
Mathematical constants
*/
namespace Constants {
    /** Names for the 9 indices of a tensor (2-rank) */
    enum {
        XX, YY, ZZ, XY, YZ, XZ, YX, ZY, ZX
    };
    const Int MAX_INT = Int(1 << 31);
    const Scalar PI = Scalar(3.14159265358979323846264);
    const Scalar E  = Scalar(2.71828182845904523536028);
#ifdef USE_DOUBLE
    const Scalar MachineEpsilon = Scalar(1e-15);
    const Scalar EqualEpsilon = Scalar(1e-7);
#else
    const Scalar MachineEpsilon = Scalar(1e-7);
    const Scalar EqualEpsilon = Scalar(1e-3);
#endif
    const Vector I_V(1,1,1);
    const Tensor I_T(1,1,1);
    const STensor I_ST(1,1,1);
#pragma acc declare create(I_V,I_T,I_ST)
}
/** Checks if two scalars are equal within machine epsilon */
FORCEINLINE  bool equal(const Scalar& p,const Scalar& q, const Scalar tol = Constants::EqualEpsilon) { 
    Scalar delta = mag(p - q);
    return delta <= tol * mag(p) ||
           delta <= tol * mag(q);
}
/** Checks if two vectors are equal within machine epsilon */
FORCEINLINE  bool equal(const Vector& p,const Vector& q, const Scalar tol = Constants::EqualEpsilon) {
    return (equal(p[0],q[0],tol) && equal(p[1],q[1],tol) && equal(p[2],q[2],tol));
}
/** \name Symmetry boundary condition*/
//@{
FORCEINLINE Scalar sym(const Scalar& p,const Vector& n) {
    return p;
}
FORCEINLINE Vector sym(const Vector& p,const Vector& n) {
    Vector en = unit(n);
    STensor A = Constants::I_ST - mul(en);
    Vector r = dot(A,p);
    Scalar magR = mag(r);
    if(equal(magR,Scalar(0)))
        return r;
    return r * (mag(p) / magR);
}
FORCEINLINE STensor sym(const STensor& p,const Vector& n) {
    Vector en = unit(n);
    STensor A = Constants::I_ST - mul(en);
    STensor r = mul(mul(A,p),A);
    Scalar magR = mag(r);
    if(equal(magR,Scalar(0)))
        return r;
    return r * (mag(p) / magR);
}
FORCEINLINE Tensor sym(const Tensor& p,const Vector& n) {
    Vector en = unit(n);
    Tensor A = Constants::I_T - mul(en,en);
    Tensor r = mul(mul(A,p),A);
    Scalar magR = mag(r);
    if(equal(magR,Scalar(0)))
        return r;
    return r * (mag(p) / magR);
}
//@}

/**
 Transfinite interpolation on face (2D)
*/
template <class T>
T Interpolate_face (Scalar r,Scalar s, T x00, T x01, T x10,
  T x11, T xr0, T xr1, T x0s, T x1s
  ) {

  T result =
     - ( 1.0 - r ) * ( 1.0 - s ) * x00
     + ( 1.0 - r )               * x0s
     - ( 1.0 - r ) *         s   * x01
     +               ( 1.0 - s ) * xr0
     +                       s   * xr1
     -         r   * ( 1.0 - s ) * x10
     +         r                 * x1s
     -         r   *         s   * x11;

  return result;
}
/**
 Transfinite interpolation on cell (3D)
*/
template <class T>
T Interpolate_cell ( Scalar r, Scalar s, Scalar t,
  T x000, T x001, T x010, T x011,
  T x100, T x101, T x110, T x111,
  T xr00, T xr01, T xr10, T xr11,
  T x0s0, T x0s1, T x1s0, T x1s1,
  T x00t, T x01t, T x10t, T x11t,
  T x0st, T x1st, T xr0t, T xr1t, T xrs0, T xrs1
  ) {

 T result =    
         ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t ) * x000
       - ( 1.0 - r ) * ( 1.0 - s )               * x00t
       + ( 1.0 - r ) * ( 1.0 - s ) *         t   * x001
       - ( 1.0 - r )               * ( 1.0 - t ) * x0s0
       + ( 1.0 - r )                             * x0st
       - ( 1.0 - r )               *         t   * x0s1
       + ( 1.0 - r ) *         s   * ( 1.0 - t ) * x010
       - ( 1.0 - r ) *         s                 * x01t
       + ( 1.0 - r ) *         s   *         t   * x011
       -               ( 1.0 - s ) * ( 1.0 - t ) * xr00
       +               ( 1.0 - s )               * xr0t
       -               ( 1.0 - s ) *         t   * xr01
       +                             ( 1.0 - t ) * xrs0
       +                                     t   * xrs1
       -                       s   * ( 1.0 - t ) * xr10
       +                       s                 * xr1t
       -                       s   *         t   * xr11
       +         r   * ( 1.0 - s ) * ( 1.0 - t ) * x100
       -         r   * ( 1.0 - s )               * x10t
       +         r   * ( 1.0 - s ) *         t   * x101
       -         r                 * ( 1.0 - t ) * x1s0
       +         r                               * x1st
       -         r                 *         t   * x1s1
       +         r   *         s   * ( 1.0 - t ) * x110
       -         r   *         s                 * x11t
       +         r   *         s   *         t   * x111;

  return result;
}

/*
 * end
 */
#endif
