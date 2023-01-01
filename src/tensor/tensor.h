#ifndef __TENSOR_H
#define __TENSOR_H

#include <iostream>
#include <cmath>

#include "my_types.h"

/**************************
 * scalars
 *************************/
#if defined __DOUBLE
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
#define Op(name,$)                                                  \
    TTensor& operator $(const TTensor& q) {                         \
        Unroll<SIZE>::name(P,q.P);                                  \
        return *this;                                               \
    }
#define SOp(name,$)                                                 \
    TTensor& operator $(const Scalar& q) {                          \
        Unroll<SIZE>::name(P,q);                                    \
        return *this;                                               \
    }
#define Fp(name)                                                    \
    friend TTensor name(const TTensor& p,const TTensor& s) {        \
        TTensor r;                                                  \
        Unroll<SIZE>::name(r.P,p.P,s.P);                            \
        return r;                                                   \
    }
#define Fp1(name)                                                   \
    friend TTensor name(const TTensor& p,const Scalar& s) {         \
        TTensor r;                                                  \
        Unroll<SIZE>::name(r.P,p.P,s);                              \
        return r;                                                   \
    }
#define Fp2(name)                                                   \
    friend TTensor name(const TTensor& p) {                         \
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
    friend Scalar dot(const TTensor& p,const TTensor& q) {
        return p & q;
    }
    friend TTensor mul(const TTensor& p,const Scalar& q) {
        return p * q;
    }
    friend Scalar magSq(const TTensor& p) {
        return p & p;
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
    friend TTensor dev(const TTensor& p,const Scalar& factor) {
        TTensor r = p;
        Scalar t = tr(p) * factor / 3;
        r[0] -= t;
        r[1] -= t;
        r[2] -= t;
        return r;
    }
    friend TTensor hyd(const TTensor& p,const Scalar& factor) {
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

/*Tensor operations*/
Vector operator ^ (const Vector& p,const Vector& q);
STensor mul(const Vector& p);
Tensor mul(const Vector& p,const Vector& q);
Tensor mul(const Tensor& p,const Tensor& q);
STensor mul(const STensor& p,const STensor& q);
Vector dot(const Tensor&,const Vector&);
Vector dot(const STensor&,const Vector&);
STensor sym(const Tensor& p);
Tensor skw(const Tensor& p);
Tensor trn(const Tensor& p);
Scalar det(const Tensor& p);
Tensor inv(const Tensor& p);
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
    const Scalar MachineEpsilon = (sizeof(Scalar) == 4) ? Scalar(1e-8) : Scalar(1e-15);
    const Scalar EqualEpsilon = (Constants::MachineEpsilon * pow(10.0,double(sizeof(Scalar))));
    const Vector I_V = Vector(1,1,1);
    const Tensor I_T = Tensor(1,1,1);
    const STensor I_ST = STensor(1,1,1);
}
/** Checks if two scalars are equal within machine epsilon */
FORCEINLINE  bool equal(const Scalar& p,const Scalar& q) { 
    return (mag(p - q) <= Constants::EqualEpsilon); 
}
/** Checks if two vectors are equal within machine epsilon */
FORCEINLINE  bool equal(const Vector& p,const Vector& q) {
    return (equal(p[0],q[0]) && equal(p[1],q[1]) && equal(p[2],q[2]));
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

/** \name Tensor array expressions */

//@{
#define DEFINE_UNARY_TOP_PART1_A(DApOpp,func,type,rtype)		                \
	class DApOpp {												                \
	public:														                \
		DApOpp() { }											                \
		static inline rtype apply(type a)						                \
		{ return func(a); }										                \
	};
#define DEFINE_UNARY_TOP_PART1_B(DApOpp,$,type,rtype)		                    \
	class DApOpp {												                \
	public:														                \
		DApOpp() { }											                \
		static inline rtype apply(type a)						                \
		{ return ($ a); }										                \
	};
#define DEFINE_UNARY_TOP_PART2(DApOpp,func,type,rtype)                          \
	DVExpr<rtype,DVUnaryExpr<rtype,DVExpr<type,A>,DApOpp> >		                \
	func(const DVExpr<type,A>& a) {						                        \
		typedef DVUnaryExpr<rtype,DVExpr<type,A>,DApOpp> ExprT;	                \
		return DVExpr<rtype,ExprT>(ExprT(a));					                \
	}

#define DEFINE_UNARY_OP(x,func,type)                                            \
    template<class type>                                                        \
    DEFINE_UNARY_TOP_PART1_A(x,func,type,type)                                  \
    template<class type, class A>								                \
    DEFINE_UNARY_TOP_PART2(x<type>,func,type,type)
        
#define DEFINE_UNARY_OP2(x,$,func,type)                                         \
    template<class type>                                                        \
    DEFINE_UNARY_TOP_PART1_B(x,$,type,type)                                     \
    template<class type, class A>								                \
    DEFINE_UNARY_TOP_PART2(x<type>,func,type,type)
        
#define DEFINE_UNARY_SOP(x,func,type,rtype)                                     \
    template<class type>                                                        \
    DEFINE_UNARY_TOP_PART1_A(x,func,type,rtype)                                 \
    template<class type, class A>                                               \
    DEFINE_UNARY_TOP_PART2(x<type>,func,type,rtype)
	

#define DEFINE_BINARY_TOP_PART1_A(DApOpp,$,func,type1,type2,rtype)		        \
	class DApOpp {												                \
	public:														                \
		DApOpp() { }											                \
		static inline rtype apply(type1 a, type2 b)				                \
		{ return (a $ b); }				                                        \
	};
#define DEFINE_BINARY_TOP_PART1_B(DApOpp,$,func,type1,type2,rtype)		        \
	class DApOpp {												                \
	public:														                \
		DApOpp() { }											                \
		static inline rtype apply(type1 a, type2 b)				                \
		{ return func(a, b); }				                                    \
	};	
#define DEFINE_BINARY_TOP_PART2(DApOpp,$,func,type1,type2,rtype)	            \
	DVExpr<rtype,DVBinExpr<rtype,DVExpr<type1,A>,DVExpr<type2,B>,DApOpp> >	    \
	func(const DVExpr<type1,A>& a, const DVExpr<type2,B>& b) {		            \
		typedef DVBinExpr<rtype,DVExpr<type1,A>,DVExpr<type2,B>,DApOpp> ExprT;	\
		return DVExpr<rtype,ExprT>(ExprT(a,b));					                \
	}
	
#define DEFINE_BINARY_OP(x,$,func,type)                                         \
    template<class type>                                                        \
    DEFINE_BINARY_TOP_PART1_A(x,$,func,type,type,type)                          \
    template<class type, class A, class B>                                      \
    DEFINE_BINARY_TOP_PART2(x<type>,$,func,type,type,type)	

#define DEFINE_BINARY_OP2(x,func,type)                                          \
    template<class type>                                                        \
    DEFINE_BINARY_TOP_PART1_B(x,$,func,type,type,type)                          \
    template<class type, class A, class B>                                      \
    DEFINE_BINARY_TOP_PART2(x<type>,$,func,type,type,type)	
   
#define DEFINE_BINARY_OP_A(x,$,func,type,rtype)                                 \
    template<class type>                                                        \
    DEFINE_BINARY_TOP_PART1_A(x,$,func,type,type,rtype)                         \
    template<class type, class A, class B>                                      \
    DEFINE_BINARY_TOP_PART2(x<type>,$,func,type,type,rtype)	

#define DEFINE_BINARY_OP2_A(x,func,type,rtype)                                  \
    template<class type>                                                        \
    DEFINE_BINARY_TOP_PART1_B(x,$,func,type,type,rtype)                         \
    template<class type, class A, class B>                                      \
    DEFINE_BINARY_TOP_PART2(x<type>,$,func,type,type,rtype)	
             
#define DEFINE_BINARY_OP_O(x,$,func,type1,type2,rtype)                          \
    DEFINE_BINARY_TOP_PART1_A(x,$,func,type1,type2,rtype)                       \
    template<class A, class B>									                \
    DEFINE_BINARY_TOP_PART2(x,$,func,type1,type2,rtype)

#define DEFINE_BINARY_OP2_O(x,func,type1,type2,rtype)                           \
    DEFINE_BINARY_TOP_PART1_B(x,$,func,type1,type2,rtype)                       \
    template<class A, class B>									                \
    DEFINE_BINARY_TOP_PART2(x,$,func,type1,type2,rtype)

#define DEFINE_BINARY_TSOP_PART1_A(DApOpp,$,func,type,type1,type2)              \
	class DApOpp {												                \
	public:														                \
		DApOpp() { }											                \
		static inline type apply(type1 a, type2 b)				                \
		{ return (a $ b); }				                                        \
	};
#define DEFINE_BINARY_TSOP_PART1_B(DApOpp,$,func,type,type1,type2)              \
	class DApOpp {												                \
	public:														                \
		DApOpp() { }											                \
		static inline type apply(type1 a, type2 b)				                \
		{ return func(a, b); }				                                    \
	};	
#define DEFINE_BINARY_TSOP_PART2_I(DApOpp,$,func,type,type2)		            \
	DVExpr<type,DVBinScaExpr<type,DVExpr<type,A>,type2,DApOpp> >	            \
	func(const DVExpr<type,A>& a, const type2& b) {			                    \
		typedef DVBinScaExpr<type,DVExpr<type,A>,type2,DApOpp> ExprT;	        \
		return DVExpr<type,ExprT>(ExprT(a,b));						            \
	}
#define DEFINE_BINARY_TSOP_PART2_II(DApOpp,$,func,type,type2)		            \
	DVExpr<type,DVBinScaInvExpr<type,type2,DVExpr<type,A>,DApOpp> >	            \
	func(const type2& a, const DVExpr<type,A>& b) {			                    \
		typedef DVBinScaInvExpr<type,type2,DVExpr<type,A>,DApOpp> ExprT;	    \
		return DVExpr<type,ExprT>(ExprT(a,b));						            \
	}
					
#define DEFINE_BINARY_SOP(x,$,func,type,type2)                                  \
    template<class type>                                                        \
    DEFINE_BINARY_TSOP_PART1_A(x,$,func,type,type,type2)                        \
    template<class type,class A>                                                \
    DEFINE_BINARY_TSOP_PART2_I(x<type>,$,func,type,type2)                       \
    template<class type>                                                        \
    DEFINE_BINARY_TSOP_PART1_A(x##Inv,$,func,type,type2,type)                   \
    template<class type,class A>                                                \
    DEFINE_BINARY_TSOP_PART2_II(x##Inv<type>,$,func,type,type2)
       
#define DEFINE_BINARY_SOP2(x,func,type,type2)                                   \
    template<class type>                                                        \
    DEFINE_BINARY_TSOP_PART1_B(x,$,func,type,type,type2)                        \
    template<class type,class A>                                                \
    DEFINE_BINARY_TSOP_PART2_I(x<type>,$,func,type,type2)                       \
    template<class type>                                                        \
    DEFINE_BINARY_TSOP_PART1_B(x##Inv,$,func,type,type2,type)                   \
    template<class type,class A>                                                \
    DEFINE_BINARY_TSOP_PART2_II(x##Inv<type>,$,func,type,type2)
        
template<class C, class A, class Op>
class DVUnaryExpr {
protected:
    A iter_;
public:
    DVUnaryExpr(const A& a)
        : iter_(a)
    { }
    C operator[](Int i) const { 
    	return Op::apply(iter_[i]); 
    }
};

template<class C, class A, class B, class Op>
class DVBinExpr {
protected:
    A iter1_;
    B iter2_;
public:
    DVBinExpr(const A& a, const B& b)
        : iter1_(a), iter2_(b)
    { }
	C operator[](Int i) const { 
    	return Op::apply(iter1_[i], iter2_[i]); 
    }
};

template<class C, class A, class B, class Op>
class DVBinScaExpr {
protected:
    A iter1_;
    B iter2_;
public:
    DVBinScaExpr(const A& a, const B& b)
        : iter1_(a), iter2_(b)
    { }
	C operator[](Int i) const { 
    	return Op::apply(iter1_[i], iter2_); 
    }
};

template<class C, class A, class B, class Op>
class DVBinScaInvExpr {
protected:
    A iter1_;
    B iter2_;
public:
    DVBinScaInvExpr(const A& a, const B& b)
        : iter1_(a), iter2_(b)
    { }
	C operator[](Int i) const { 
    	return Op::apply(iter1_, iter2_[i]); 
    }
};

template<class type, class C>
class DVExpr {
protected:
    C P;
public:
    DVExpr()
    { }
	DVExpr(const C& a)
		: P(a)
	{ }
    type operator[](Int i) const
    { return P[i]; }
};

//templated expressions
DEFINE_UNARY_OP(oacos,acos,type);
DEFINE_UNARY_OP(oasin,asin,type);
DEFINE_UNARY_OP(oatan,atan,type);
DEFINE_UNARY_OP(oceil,ceil,type);
DEFINE_UNARY_OP(ocos,cos,type);
DEFINE_UNARY_OP(ocosh,cosh,type);
DEFINE_UNARY_OP(oexp,exp,type);
DEFINE_UNARY_OP(ofabs,fabs,type);
DEFINE_UNARY_OP(ofloor,floor,type);
DEFINE_UNARY_OP(olog,log,type);
DEFINE_UNARY_OP(olog10,log10,type);
DEFINE_UNARY_OP(osin,sin,type);
DEFINE_UNARY_OP(osinh,sinh,type);
DEFINE_UNARY_OP(osqrt,sqrt,type);
DEFINE_UNARY_OP(otan,tan,type);
DEFINE_UNARY_OP(otanh,tanh,type);

DEFINE_UNARY_OP(otrn,trn,Tensor);
DEFINE_UNARY_OP(oskw,skw,Tensor);
DEFINE_UNARY_OP(ounit,unit,type);

DEFINE_UNARY_OP2(oneg,-,operator -,type);

DEFINE_UNARY_SOP(omag,mag,type,Scalar);
DEFINE_UNARY_SOP(osym,sym,Tensor,STensor)
    
DEFINE_BINARY_OP(oadd,+,operator +,type);
DEFINE_BINARY_OP(osub,-,operator -,type);
DEFINE_BINARY_OP(omul,*,operator *,type);
DEFINE_BINARY_OP(odiv,/,operator /,type);

DEFINE_BINARY_OP2(oatan2,atan2,type);
DEFINE_BINARY_OP2(omin,min,type);
DEFINE_BINARY_OP2(omax,max,type);
DEFINE_BINARY_OP2(ossdiv,sdiv,type);

DEFINE_BINARY_SOP(osadd,+,operator +,type,Scalar);
DEFINE_BINARY_SOP(ossub,-,operator -,type,Scalar);
DEFINE_BINARY_SOP(osmul,*,operator *,type,Scalar);
DEFINE_BINARY_SOP(osdiv,/,operator /,type,Scalar);

DEFINE_BINARY_SOP2(opow,pow,type,Scalar);
DEFINE_BINARY_SOP2(odev,dev,type,Scalar);
DEFINE_BINARY_SOP2(ohyd,hyd,type,Scalar);
DEFINE_BINARY_SOP2(omins,min,type,type);
DEFINE_BINARY_SOP2(omaxs,max,type,type);

DEFINE_BINARY_OP_A(oinner,&,operator &,type,Scalar)
DEFINE_BINARY_OP2_A(odot,dot,type,Scalar)
    
//specific expressions
#define ScaMul(T,$)                                          \
    DEFINE_BINARY_OP_O(omul1##$,*,operator *,T,Scalar,T)     \
    DEFINE_BINARY_OP_O(omul2##$,*,operator *,Scalar,T,T)     \
    DEFINE_BINARY_OP_O(odiv1##$,/,operator /,T,Scalar,T)     \
    DEFINE_BINARY_OP_O(odiv2##$,/,operator /,Scalar,T,T)

    ScaMul(Vector,1)
    ScaMul(Tensor,2)
    ScaMul(STensor,3)

#undef ScaMul

DEFINE_BINARY_OP_O(omula,*,operator *,Scalar,Scalar,Scalar);
DEFINE_BINARY_OP_O(omulb,*,operator *,Vector,Vector,Vector);
DEFINE_BINARY_OP_O(omulc,*,operator *,Tensor,Tensor,Tensor);
DEFINE_BINARY_OP_O(omuld,*,operator *,STensor,STensor,STensor);

DEFINE_BINARY_OP2_O(omul2,mul,Vector,Vector,Tensor)
DEFINE_BINARY_OP2_O(omul3,mul,Vector,Scalar,Vector)
DEFINE_BINARY_OP2_O(omul4,mul,Tensor,Tensor,Tensor)
DEFINE_BINARY_OP2_O(omul5,mul,STensor,STensor,STensor)

DEFINE_BINARY_OP2_O(odot2,dot,STensor,Vector,Vector)
DEFINE_BINARY_OP2_O(odot3,dot,Tensor,Vector,Vector)
//@}
    
/*
 * end
 */
#endif
