#include "solve.h"

/* *********************************************************************
 *  Solve system of linear equations iteratively
 * *********************************************************************/
template<class type>
Scalar getResidual(const MeshField<type,CELL>& r,
				   const MeshField<type,CELL>& cF,
				   bool sync) {
	type res[2];
	res[0] = type(0);
	res[1] = type(0); 
	for(Int i = 0;i < Mesh::gBCellsStart;i++) {
		res[0] += (r[i] * r[i]);
		res[1] += (cF[i] * cF[i]);
	}
	if(sync) {
		type global_res[2];
		MP::allsum(res,global_res,2);
		res[0] = global_res[0];
		res[1] = global_res[1];
	}
	return sqrt(sdiv(mag(res[0]), mag(res[1])));
}

template<class type>
void SolveT(const MeshMatrix<type>& M) {

	using namespace Mesh;
	MeshField<type,CELL> r,p,AP;
	MeshField<type,CELL> r1(false),p1(false),AP1(false);   
	MeshField<type,CELL>& cF = *M.cF;
	MeshField<type,CELL>& buffer = AP;
	ScalarCellField D = M.ap,iD = (1 / M.ap);
	Scalar res,ires;
	type alpha,beta,o_rr = type(0),oo_rr;
	Int iterations = 0;
	bool converged = false;
	register Int i;

	/****************************
	 * Parallel controls
	 ***************************/
	bool print = (MP::host_id == 0);
	int  end_count = 0;
	bool sync = (Controls::parallel_method == Controls::BLOCKED)
		&& gInterMesh.size();
	std::vector<bool> sent_end(gInterMesh.size(),false);

	/****************************
	 * Identify solver type
	 ***************************/
	if(print) {
		if(M.flags & M.SYMMETRIC)
			MP::printH("SYMM-");
		else
			MP::printH("ASYM-");
		if(Controls::Solver == Controls::JACOBI)
			MP::print("JAC :");
		else if(Controls::Solver == Controls::SOR)
			MP::print("SOR :");
		else {
			switch(Controls::Preconditioner) {
			case Controls::NOP: MP::print("PCG :"); break;
			case Controls::DIAG: MP::print("DIAG-PCG :"); break;
			case Controls::SORP: MP::print("SOR-PCG :"); break;
			case Controls::DILU: MP::print("DILU-PCG :"); break;
			}
		}
	}
	/****************************
	 * Initialization
	 ***************************/
	if(Controls::Solver == Controls::PCG) {
		if(!(M.flags & M.SYMMETRIC)) {
			/* Allocate BiCG vars*/
			r1.allocate();
			p1.allocate();
			AP1.allocate();
		} else {
			if(Controls::Preconditioner == Controls::SORP) {
				/*SOR and GS*/
				iD *= Controls::SOR_omega;
				D *=  (2.0 / Controls::SOR_omega - 1.0);	
			} else if(Controls::Preconditioner == Controls::DILU) {
				/*D-ILU(0)*/
				for(i = 0;i < gBCellsStart;i++) {
					Cell& c = gCells[i];
					forEach(c,j) {								
						Int f = c[j];							
						Int c1 = gFO[f];						
						Int c2 = gFN[f];						
						if(i == c1) {
							if(c2 > i) D[c2] -= 
								(M.an[0][f] * M.an[1][f] * iD[c1]);	
						} else {
							if(c1 > i) D[c1] -= 
								(M.an[0][f] * M.an[1][f] * iD[c2]);		
						}										
					}			
				}
				iD = (1 / D);
			}
			/*end*/
		}
	}
	/****************************
	 * Jacobi sweep
	 ***************************/
#define JacobiSweep() {								\
	AP = iD * getRHS(M);							\
	for(i = 0;i < gBCellsStart;i++)					\
		cF[i] = AP[i];								\
}
	/****************************
	 *  Forward/backward GS sweeps
	 ****************************/
#define Sweep_(X,B,i) {								\
	Cell& c = gCells[i];							\
	type ncF = B[i];								\
	forEach(c,j) {									\
		Int f = c[j];								\
		if(i == gFO[f])								\
			ncF += X[gFN[f]] * M.an[1][f];			\
		else										\
			ncF += X[gFO[f]] * M.an[0][f];			\
	}												\
	ncF *= iD[i];									\
	X[i] = X[i] * (1 - Controls::SOR_omega) +		\
		ncF * (Controls::SOR_omega);				\
}
#define ForwardSweep(X,B) {							\
	for(i = 0;i < gBCellsStart;i++)					\
		Sweep_(X,B,i);								\
}
#define BackwardSweep(X,B) {						\
	for(int i = gBCellsStart - 1;i >= 0;i--)		\
		Sweep_(X,B,i);								\
}
	/***********************************
	 *  Forward/backward substitution
	 ***********************************/
#define Substitute_(X,B,i,forw,tr) {				\
	Cell& c = gCells[i];							\
	type ncF = B[i];								\
	forEach(c,j) {									\
		Int f = c[j];								\
		Int c1 = gFO[f];							\
	    Int c2 = gFN[f];							\
		if(i == c1) {								\
			if((forw && (c2 < c1)) ||				\
			  (!forw && (c1 < c2)))	{				\
				ncF += X[c2] * M.an[1 - tr][f];		\
			}										\
		} else {									\
			if((forw && (c2 > c1)) ||				\
			  (!forw && (c1 > c2)))					\
				ncF += X[c1] * M.an[0 + tr][f];		\
		}											\
	}												\
	ncF *= iD[i];									\
	X[i] = ncF;										\
}
#define ForwardSub(X,B,TR) {						\
	for(i = 0;i < gBCellsStart;i++)					\
		Substitute_(X,B,i,true,TR);					\
}
#define BackwardSub(X,B,TR) {						\
	for(i = gBCellsStart;i > 0;i--)					\
		Substitute_(X,B,i-1,false,TR);				\
}
#define DiagSub(X,B) {								\
	for(i = 0;i < gBCellsStart;i++)					\
		X[i] = B[i] * iD[i];						\
}
	/***********************************
	 *  Preconditioners
	 ***********************************/
#define precondition_(R,Z,TR) {						\
	using namespace Controls;						\
	if(Preconditioner == Controls::NOP) {			\
		Z = R;										\
	} else if(Preconditioner == Controls::DIAG) {	\
		DiagSub(Z,R);								\
	} else {										\
		if(Controls::Solver == Controls::PCG) {		\
			Z = type(0);							\
			ForwardSub(Z,R,TR);						\
			Z = Z * D;								\
			BackwardSub(Z,Z,TR);					\
		}											\
	}												\
}
#define precondition(R,Z) precondition_(R,Z,0)
#define preconditionT(R,Z) precondition_(R,Z,1)
	/***********************************
	 *  SAXPY and DOT operations
	 ***********************************/
#define Taxpy(Y,I,X,alpha_) {						\
	for(i = 0;i < gBCellsStart;i++)					\
		Y[i] = I[i] + X[i] * alpha_;				\
}
#define Tdot(X,Y,sum) {								\
	sum = type(0);									\
	for(i = 0;i < gBCellsStart;i++)					\
		sum += X[i] * Y[i];							\
}
	/***********************************
	 *  Synchronized sum and exchange
	 ***********************************/
#define SUM_ALL(typ,var)	if(sync)	{			\
	typ t;											\
	MP::allsum(&var,&t,1);							\
	var = t;										\
}
#define EXCHANGE(var)		if(sync)	{			\
	exchange_ghost(&var[0]);						\
}
	/***********************************
	 *  Residual
	 ***********************************/
#define CALC_RESID() {								\
	r = M.Su - mul(M,cF);							\
	forEachS(r,k,gBCellsStart)						\
		r[k] = type(0);								\
	precondition(r,AP);								\
	forEachS(AP,k,gBCellsStart)						\
		AP[k] = type(0);							\
	res = getResidual(AP,cF,sync);					\
	if(Controls::Solver == Controls::PCG) {			\
		Tdot(r,AP,o_rr);							\
		SUM_ALL(type,o_rr);							\
		p = AP;										\
		if(!(M.flags & M.SYMMETRIC)) {				\
			r1 = r;									\
			p1 = p;									\
		}											\
	}												\
}
	/***********************
	 *  Initialize residual
	 ***********************/
	CALC_RESID();
	ires = res;
	/********************************************************
	* Initialize exchange of ghost cells just once.
	* Lower numbered processors send message to higher ones.
	*********************************************************/
	if(!sync) {
		end_count = gInterMesh.size();
		forEach(gInterMesh,i) {
			interBoundary& b = gInterMesh[i];
			if(b.from < b.to) {
				IntVector& f = *(b.f);
				forEach(f,j)
					buffer[j] = cF[gFO[f[j]]];
				MP::send(&buffer[0],f.size(),b.to,MP::FIELD);
			}
		}
	}
    /* **************************
	 * Iterative solution
	 * *************************/
	while(iterations < Controls::max_iterations) {
		/*counter*/
		iterations++;

		/*select solver*/
		if(Controls::Solver == Controls::JACOBI) {
			/*Jacobi solver*/
			p = cF;
			JacobiSweep();
			for(i = 0;i < gBCellsStart;i++)
				AP[i] = cF[i] - p[i];
			/*end*/
		} else if(Controls::Solver == Controls::SOR) {
			/*Asynchronous SOR solver*/
			p = cF;
			ForwardSweep(cF,M.Su);
			for(i = 0;i < gBCellsStart;i++)
				AP[i] = cF[i] - p[i];
			/*end*/
		} else if(M.flags & M.SYMMETRIC) {
			/*conjugate gradient*/
			EXCHANGE(p);
			AP = mul(M,p);
			Tdot(p,AP,oo_rr);
			SUM_ALL(type,oo_rr);
			alpha = sdiv(o_rr , oo_rr);
			Taxpy(cF,cF,p,alpha);
			Taxpy(r,r,AP,-alpha);
			precondition(r,AP);
			oo_rr = o_rr;
			Tdot(r,AP,o_rr);
			SUM_ALL(type,o_rr);
			beta = sdiv(o_rr , oo_rr);
			Taxpy(p,AP,p,beta);
			/*end*/
		} else {
			/* biconjugate gradient*/
			EXCHANGE(p);
			EXCHANGE(p1);
			AP = mul(M,p);
			AP1 = mult(M,p1);
			Tdot(p1,AP,oo_rr);
			SUM_ALL(type,oo_rr);
			alpha = sdiv(o_rr , oo_rr);
			Taxpy(cF,cF,p,alpha);
			Taxpy(r,r,AP,-alpha);
			Taxpy(r1,r1,AP1,-alpha);
			precondition(r,AP);
			preconditionT(r1,AP1);
			oo_rr = o_rr;
			Tdot(r1,AP,o_rr);
			SUM_ALL(type,o_rr);
			beta = sdiv(o_rr , oo_rr);
			Taxpy(p,AP,p,beta);
			Taxpy(p1,AP1,p1,beta);
			/*end*/
		}
		/* *********************************************
		* calculate norm of residual & check convergence
		* **********************************************/
		EXCHANGE(cF);
		res = getResidual(AP,cF,sync);
		if(res <= Controls::tolerance
			|| iterations == Controls::max_iterations)
			converged = true;
PROBE:
		/* **********************************************************
		 * Update ghost cell values. Communication is NOT forced on 
		 * every iteration,rather a non-blocking probe is used to 
		 * process messages as they arrive.
		 ************************************************************/
		if(!sync)
		{
			int source,message_id;
			/*probe*/
			while(MP::iprobe(source,message_id)) {
				/*find the boundary*/
				Int patchi;
				for(patchi = 0;patchi < gInterMesh.size();patchi++) {
					if(gInterMesh[patchi].to == (Int)source) 
						break;
				}
				interBoundary& b = gInterMesh[patchi];
				/*parse message*/
				if(message_id == MP::FIELD) {
					IntVector& f = *(b.f);
					/*recieve*/
					MP::recieve(&buffer[0],f.size(),source,message_id);
					forEach(f,j)
						cF[gFN[f[j]]] = buffer[j];
					/*Re-calculate residual.*/
					CALC_RESID();
					if(res > Controls::tolerance
						&& iterations < Controls::max_iterations)
						converged = false;
					/* For communication to continue, processor have to send back 
					 * something for every message recieved.*/
					if(converged) {
						/*send END marker*/
						if(!sent_end[patchi]) {
							MP::send(source,MP::END);
							sent_end[patchi] = true;
						}
					} else {
						/*send back our part*/
						forEach(f,j)
							buffer[j] = cF[gFO[f[j]]];
						MP::send(&buffer[0],f.size(),source,message_id);
					}
				} else if(message_id == MP::END) {
					/*END marker recieved*/
					MP::recieve(source,message_id);
					end_count--;
					if(!sent_end[patchi]) {
						MP::send(source,MP::END);
						sent_end[patchi] = true;
					}
				}
			}
		}
		/* *****************************************
		* Wait untill all partner processors send us
		* an END message i.e. until end_count = 0.
		* *****************************************/
		if(converged) {
			if(end_count > 0) goto PROBE;
			else break;
		}
		/********
		 * end
		 ********/
	}

	/*solver info*/
	if(print)
		MP::print("Iterations %d Initial Residual "
		"%.5e Final Residual %.5e\n",iterations,ires,res);

	/*barrier*/
	MP::barrier();

	/*update boundary conditons*/
	updateExplicitBCs(cF);
}
/***************************
 * Explicit instantiations
 ***************************/
void Solve(const MeshMatrix<Scalar>& A) {
	applyImplicitBCs(A);
	SolveT(A);
}
void Solve(const MeshMatrix<Vector>& A) {
	applyImplicitBCs(A);
	SolveT(A);
}
void Solve(const MeshMatrix<STensor>& A) {
	applyImplicitBCs(A);
	SolveT(A);
}
void Solve(const MeshMatrix<Tensor>& A) {
	applyImplicitBCs(A);
	SolveT(A);
}
/* ********************
 *        End
 * ********************/
