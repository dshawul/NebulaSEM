#ifndef __DG_H
#define __DG_H

#define forEachLgl(i_,j_,k_) 							\
			for(Int i_ = 0;i_ < NPX;i_++)				\
				for(Int j_ = 0;j_ < NPY;j_++)			\
					for(Int k_ = 0;k_ < NPZ;k_++)
#define forEachLglXY(i_,j_) 							\
			for(Int i_ = 0;i_ < NPX;i_++)				\
				for(Int j_ = 0;j_ < NPY;j_++)
#define forEachLglXZ(i_,k_) 							\
			for(Int i_ = 0;i_ < NPX;i_++)				\
				for(Int k_ = 0;k_ < NPZ;k_++)
#define forEachLglYZ(j_,k_) 							\
			for(Int j_ = 0;j_ < NPY;j_++)				\
				for(Int k_ = 0;k_ < NPZ;k_++)					
				
#define INDEX4(c_,i_,j_,k_) \
	((c_) * NPX * NPY * NPZ + (i_) * NPY * NPZ + (j_) * NPZ + (k_))
	
#define INDEX3(i_,j_,k_) \
	((i_) * NPY * NPZ + (j_) * NPZ + (k_))
	
namespace DG {
	extern Scalar **psi[3];
	extern Scalar **dpsi[3];
	extern Scalar *xgl[3];
	extern Scalar *wgl[3];
	extern TensorCellField Jinv;
	extern Int NPX,NPY,NPZ;
	
	void legendre(int p, Scalar x,Scalar& L0,Scalar& L0_1,Scalar& L0_2);
	void legendre_gauss_lobatto(int N, Scalar* xgl, Scalar* wgl);
	void lagrange(Scalar x,int N, Scalar* xgl,Scalar* psi);
	void lagrange_der(Scalar x,int N, Scalar* xgl,Scalar* dpsi);
	void init_poly();
	void init_basis();
	void init_geom();
	
	/*expand fields*/
	template<class type, ENTITY entity>
	void expand(const MeshField<type,entity>& cF) {
		if(DG::NPMAT) {
			Int block;
			if(entity == CELL) 
				block = NP;
			else if(entity == FACET) 
				block = NPF;
			
			for(int i = cF.size() - 1;i >= 0;i -= block) {
				Int ii = i / block;
				type C = cF[ii];
				for(Int j = 0; j < block;j++)
					cF[i - j] = C; 
			}
		}
	}
}

#endif