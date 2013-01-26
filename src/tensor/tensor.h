#ifndef __TENSOR_H
#define __TENSOR_H

#ifdef _MSC_VER
#    pragma warning (disable: 4996)
#endif

#include <fstream>
#include <iostream>
#include <cmath>

#define __DOUBLE

#ifdef _MSC_VER
#	define FORCEINLINE __forceinline
#else
#	define FORCEINLINE __inline
#endif

/*******************
 * Int is unsigned
 ******************/
typedef unsigned   Int;

/**************************
 * scalars
 *************************/
#if defined __DOUBLE
#	define Scalar  double
#else
#	define Scalar  float
#endif

/* Arthimetic operators are defined via compound assignment*/
#define Operator(T,$)												\
	friend FORCEINLINE T operator $ (const T& p,const T& q) {		\
		T r = p;													\
		r $##= q;													\
		return r;													\
	}
/* Default operator overloads for scalars vs others*/
#define OpS($)														\
	template<class T>												\
	FORCEINLINE T operator $ (const T& p,const Scalar& q) {			\
		T r = p;													\
		r $##= q;													\
		return r;													\
	}
#define COp($)														\
	template<class T>												\
	FORCEINLINE T operator $ (const Scalar& p,const T& q) {			\
		T r = q;													\
		r $##= p;													\
		return r;													\
	}
#define NCOp($)														\
	template<class T>												\
	FORCEINLINE T operator $ (const Scalar& p,const T& q) {			\
		T r = p;													\
		r $##= q;													\
		return r;													\
	}
OpS(+);
OpS(-);
OpS(*);
OpS(/);
COp(+);
COp(*);
NCOp(/);
NCOp(-);
#undef OpS
#undef COp
#undef NCOp
/*other scalar operations*/
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
/*********************************
 * loop unroller for tensors
 *********************************/
template <int N>
struct Unroll {
	/*macro*/
#define Op(name,$)																\
	static FORCEINLINE void name(Scalar* p,const Scalar* q) {					\
		*p $ *q;																\
		Unroll<N - 1>::name(p + 1,q + 1);										\
	}     
#define SOp(name,$)																\
	static FORCEINLINE void name(Scalar* p,const Scalar q) {					\
		*p $ q;																	\
		Unroll<N - 1>::name(p + 1,q);											\
	}
#define Fp(name,$)                                                              \
	static FORCEINLINE void name(Scalar* r,const Scalar* p,const Scalar* q) {	\
		*r = ::$(*p,*q);												        \
		Unroll<N - 1>::name(r + 1,p + 1,q + 1);					                \
	}
#define Fp1(name,$)                                                             \
	static FORCEINLINE void name(Scalar* r,const Scalar* p,const Scalar q) {	\
		*r = ::$(*p,q);															\
		Unroll<N - 1>::name(r + 1,p + 1,q);										\
	}
#define Fp2(name,$)                                                             \
	static FORCEINLINE void name(Scalar* r,const Scalar* p) {					\
		*r = ::$(*p);													        \
		Unroll<N - 1>::name(r + 1,p + 1);						                \
	}
	/*special*/
	static FORCEINLINE Scalar dot(const Scalar* p,const Scalar* q) {	 
		return (*p) * (*q) + Unroll<N - 1>::dot(p + 1,q + 1);					                 
	}
	/*define ops*/
	Op(equ,=);
	Op(neg,=-);
	Op(inc,+=);
	Op(dec,-=);
	Op(mul,*=);
	Op(div,/=);
	SOp(equ,=);
	SOp(neg,=-);
	SOp(inc,+=);
	SOp(dec,-=);
	SOp(mul,*=);
	SOp(div,/=);
	Fp(sdiv,sdiv);
    /*from math.h*/
	Fp2(acos,acos);
	Fp2(asin,asin);
	Fp2(atan,atan);
	Fp(atan2,atan2);
	Fp2(ceil,ceil);
	Fp2(cos,cos);
	Fp2(cosh,cosh);
	Fp2(exp,exp);
	Fp2(fabs,fabs);
	Fp2(floor,floor);
	Fp2(log,log);
	Fp2(log10,log10);
	Fp1(pow,pow);
	Fp2(sin,sin);
	Fp2(sinh,sinh);
	Fp2(sqrt,sqrt);
	Fp2(tan,tan);
	Fp2(tanh,tanh);
	Fp(min,min);
	Fp(max,max);
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
};

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
	Op(equ);
	Op(neg);
	Op(inc);
	Op(dec);
	Op(mul);
	Op(div);
	SOp(equ);
	SOp(neg);
	SOp(inc);
	SOp(dec);
	SOp(mul);
	SOp(div);
	Fp(sdiv);
	/*from math.h*/
	Fp2(acos);
	Fp2(asin);
	Fp2(atan);
	Fp(atan2);
	Fp2(ceil);
	Fp2(cos);
	Fp2(cosh);
	Fp2(exp);
	Fp2(fabs);
	Fp2(floor);
	Fp2(log);
	Fp2(log10);
	Fp1(pow);
	Fp2(sin);
	Fp2(sinh);
	Fp2(sqrt);
	Fp2(tan);
	Fp2(tanh);
	Fp(min);
	Fp(max);
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
};

/***************************************
 * Template Tensor class
 ***************************************/
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
	Scalar& operator [] (Int i) {
		return P[i];
	}
	const Scalar& operator [] (Int i) const {
		return P[i];
	}
	/*unary ops*/
	TTensor operator - () {
		TTensor r;
		Unroll<SIZE>::neg(r.P,P);
		return r;
	}
	friend Scalar operator & (const TTensor& p,const TTensor& q) {
		Scalar r = Unroll<SIZE>::dot(p.P,q.P);
		if(SIZE == 6) r += Unroll<3>::dot(&p.P[3],&q.P[3]);
		return r;
	}

	/*unrolled operations*/
#define Op(name,$)													\
	TTensor& operator $(const TTensor& q) {							\
		Unroll<SIZE>::name(P,q.P);									\
		return *this;												\
	}
#define SOp(name,$)													\
	TTensor& operator $(const Scalar& q) {							\
		Unroll<SIZE>::name(P,q);									\
		return *this;												\
	}
#define Fp(name)													\
	friend TTensor name(const TTensor& p,const TTensor& s) {		\
		TTensor r;													\
		Unroll<SIZE>::name(r.P,p.P,s.P);							\
		return r;													\
	}
#define Fp1(name)													\
	friend TTensor name(const TTensor& p,const Scalar& s) {			\
		TTensor r;													\
		Unroll<SIZE>::name(r.P,p.P,s);								\
		return r;													\
	}
#define Fp2(name)													\
	friend TTensor name(const TTensor& p) {							\
		TTensor r;													\
		Unroll<SIZE>::name(r.P,p.P);								\
		return r;													\
	}
    /*define ops*/
	Op(equ,=);
	Op(inc,+=);
	Op(dec,-=);
	Op(mul,*=);
	Op(div,/=);
	SOp(equ,=);
	SOp(inc,+=);
	SOp(dec,-=);
	SOp(mul,*=);
	SOp(div,/=);
	Fp(sdiv);
	/*from math.h*/
	Fp2(acos);
	Fp2(asin);
	Fp2(atan);
	Fp(atan2);
	Fp2(ceil);
	Fp2(cos);
	Fp2(cosh);
	Fp2(exp);
	Fp2(fabs);
	Fp2(floor);
	Fp2(log);
	Fp2(log10);
	Fp1(pow);
	Fp2(sin);
	Fp2(sinh);
	Fp2(sqrt);
	Fp2(tan);
	Fp2(tanh);
	Fp(min);
	Fp(max);
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
	Operator(TTensor,+);
	Operator(TTensor,-);
	Operator(TTensor,*);
	Operator(TTensor,/);
	/*others*/
	friend Scalar magSq(const TTensor& p) {
		return (p & p);
	}
	friend Scalar mag(const TTensor& p) {
		return sqrt(magSq(p));
	}
	friend TTensor unit(const TTensor& p) {
		TTensor r = p;
		Scalar mg = mag(r);
		r /= mg;
		return r;
	}
	friend Scalar tr(const TTensor& p) {
		return p[0] + p[1] + p[2];
	}
	friend TTensor dev(const TTensor& p,const Scalar factor = 1.) {
		TTensor r = p;
		Scalar t = tr(p) * factor / 3;
		r[0] -= t;
		r[1] -= t;
		r[2] -= t;
		return r;
	}
	friend TTensor hyd(const TTensor& p,const Scalar factor = 1.) {
		TTensor r(1,1,1);
		Scalar t = tr(p) * factor / 3;
		r[0] = t;
		r[1] = t;
		r[2] = t;
		return r;
	}
	/*IO*/
	friend std::ostream& operator << (std::ostream& os, const TTensor<SIZE>& p) {
		for(Int i = 0;i < SIZE;i++) 
			os << p[i] << " ";
		return os;
	}
	friend std::istream& operator >> (std::istream& is, TTensor<SIZE>& p) {
		for(Int i = 0;i < SIZE;i++) 
			is >> p[i];
		return is;
	}
};
/*typedef tensors*/
typedef TTensor<3> Vector;
typedef TTensor<6> STensor;
typedef TTensor<9> Tensor;

/*Tensor operations*/
Vector operator ^ (const Vector& p,const Vector& q);
STensor mul(const Vector& p);
Tensor mul(const Vector& p,const Vector& q);
Tensor mul(const Tensor& p,const Tensor& q);
STensor mul(const STensor& p,const STensor& q);
Vector dot(const Vector&,const Tensor&);
Vector dot(const Vector&,const STensor&);
STensor sym(const Tensor& p);
Tensor skw(const Tensor& p);
Tensor trn(const Tensor& p);

/*constants*/
namespace Constants {
	enum {
		XX, YY , ZZ , XY, YZ , XZ , YX, ZY , ZX
	};
	const Int MAX_INT = Int(1 << 31);
	const Scalar PI = Scalar(3.14159265358979323846264);
	const Scalar E  = Scalar(2.71828182845904523536028);
	const Scalar MachineEpsilon = (sizeof(Scalar) == 4) ? Scalar(1e-8) : Scalar(1e-15);
	const Vector I_V = Vector(1,1,1);
	const Tensor I_T = Tensor(1,1,1);
	const STensor I_ST = STensor(1,1,1);
}

FORCEINLINE  bool equal(const Scalar& p,const Scalar& q) { 
	return mag(p - q) <= (Constants::MachineEpsilon * pow(10.0,double(sizeof(Scalar)))); 
}
FORCEINLINE  bool equal(const Vector& p,const Vector& q) {
	return (equal(p[0],q[0]) && equal(p[1],q[1]) && equal(p[2],q[2]));
}
/*for symmetry boundary condition*/
FORCEINLINE Scalar sym(const Scalar& p,const Vector& n) {
	return p;
}
FORCEINLINE Vector sym(const Vector& p,const Vector& n) {
	Vector en = unit(n);
	STensor A = Constants::I_ST - mul(en);
	Vector r = dot(p,A);
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
/*
 * Blending
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
