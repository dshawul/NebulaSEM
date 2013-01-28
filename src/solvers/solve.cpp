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
	return sqrt(mag(res[0]) / mag(res[1]));
}

template<class type>
void SolveT(const MeshMatrix<type>& M) {

	using namespace Mesh;
	MeshField<type,CELL> r,p,AP;
	MeshField<type,CELL> r1(false),p1(false),AP1(false);   
	MeshField<type,CELL>& cF = *M.cF;
	MeshField<type,CELL>& buffer = AP;
	ScalarCellField pC = (1 / M.ap); /* Jacobi preconditioner */
	Scalar res,ires;
	Int i,j,iterations = 0;
	type alpha,beta,o_rr = type(0),oo_rr;

	/* Allocate BiCG vars*/
	if((Controls::Solver == Controls::PCG) &&
		!(M.flags & M.SYMMETRIC)) {
		r1.allocate();
		p1.allocate();
		AP1.allocate();
	}

	/*for cluster use*/
	bool converged = false;
	bool print = (MP::host_id == 0);
	int  end_count = 0;
	bool sync = (Controls::parallel_method == Controls::BLOCKED)
		&& gInterMesh.size();
	std::vector<bool> sent_end(gInterMesh.size(),false);

	/*identify matrix solver type*/
	if(print) {
		if(M.flags & M.SYMMETRIC)
			MP::printH("Symmetric  : ");
		else
			MP::printH("Asymmetric : ");
		if(Controls::Solver == Controls::JACOBI)
			MP::print("JAC :");
		else if(Controls::Solver == Controls::SOR)
			MP::print("SOR :");
		else
			MP::print("PCG :");
	}
	/*****************************************
	 *  Initialize residual and other vectors
	 *****************************************/
#define SUM_ALL(typ,var)	if(sync)	{		\
	typ t;										\
	MP::allsum(&var,&t,1);						\
	var = t;									\
};

#define EXCHANGE(var)		if(sync)	{		\
	exchange_ghost(&var[0]);					\
};

#define CALC_RESID() {							\
	r = M.Su - M * cF;							\
	forEachS(r,i,gBCellsStart)					\
		r[i] = type(0);							\
	res = getResidual(r * pC,cF,sync);			\
	if(Controls::Solver == Controls::PCG) {		\
		forEachS(r,i,gBCellsStart)				\
			p[i] = type(0);						\
		o_rr = type(0);							\
		for(Int i = 0;i < gBCellsStart;i++) {	\
			p[i] = r[i] * pC[i];				\
			o_rr += r[i] * p[i];				\
		}										\
		SUM_ALL(type,o_rr);						\
		if(!(M.flags & M.SYMMETRIC)) {			\
			r1 = r;								\
			p1 = p;								\
		}										\
	}											\
};

	CALC_RESID();
	ires = res;
	/************************************************
	* Initialize exchange of ghost cells just once.
	* Lower numbered processors send message to higher ones.
	************************************************/
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
			/*Jacobi solver: good for debugging*/
			p = pC * getRHS(M);
			for(i = 0;i < gBCellsStart;i++) {
				r[i] = (p[i] - cF[i]) / pC[i];
				cF[i] = p[i];
			}
		} else if(Controls::Solver == Controls::SOR) {
			/*Asynchronous SOR solver*/
			Cell* c;
			Int sz,f;
			type ncF;
			Scalar pc;
			for(i = 0;i < gBCellsStart;i++) {
				c = &gCells[i];
				sz = c->size();
				pc = pC[i];
				ncF = (M.Su[i] * pc);
				for(j = 0;j < sz;j++) {
					f = (*c)[j];
					if(i == gFO[f]) {
						ncF += cF[gFN[f]] * (M.an[1][f] * pc);
					} else {
						ncF += cF[gFO[f]] * (M.an[0][f] * pc);
					}
				}
				ncF = cF[i] * (1 - Controls::SOR_omega) + 
					  ncF * (Controls::SOR_omega);
				r[i] = (ncF - cF[i]) / pC[i];
				cF[i] = ncF;
			}
			/*end*/
		} else if(M.flags & M.SYMMETRIC) {
			/*conjugate gradient*/
			EXCHANGE(p);
			AP = M * p;
			oo_rr = type(0);
			for(i = 0;i < gBCellsStart;i++)
				oo_rr += p[i] * AP[i];
			SUM_ALL(type,oo_rr);
			alpha = sdiv(o_rr , oo_rr);
			for(i = 0;i < gBCellsStart;i++) {
				cF[i] = cF[i] + p[i] * alpha;
				r[i] = r[i] - AP[i] * alpha;
			}
			oo_rr = o_rr;
			o_rr = type(0);
			for(i = 0;i < gBCellsStart;i++) {
				AP[i] = r[i] * pC[i];
				o_rr += r[i] * AP[i];
			}
			SUM_ALL(type,o_rr);
			beta = sdiv(o_rr , oo_rr);
			for(i = 0;i < gBCellsStart;i++) {
				p[i] = AP[i] + p[i] * beta;
			}
			/*end*/
		} else {
			/* biconjugate gradient*/
			EXCHANGE(p);
			EXCHANGE(p1);
			AP = M * p;
			AP1 = M ^ p1;
			oo_rr = type(0);
			for(i = 0;i < gBCellsStart;i++)
				oo_rr += p1[i] * AP[i];
			SUM_ALL(type,oo_rr);
			alpha = sdiv(o_rr , oo_rr);
			for(i = 0;i < gBCellsStart;i++) {
				cF[i] = cF[i] + p[i] * alpha;
				r[i] = r[i] - AP[i] * alpha;
				r1[i] = r1[i] - AP1[i] * alpha;
			}
			oo_rr = o_rr;
			o_rr = type(0);
			for(i = 0;i < gBCellsStart;i++) {
				AP[i] = r[i] * pC[i];
				AP1[i] = r1[i] * pC[i];
				o_rr += r1[i] * AP[i];
			}
			SUM_ALL(type,o_rr);
			beta = sdiv(o_rr , oo_rr);
			for(i = 0;i < gBCellsStart;i++) {
				p[i] = AP[i] + p[i] * beta;
				p1[i] = AP1[i] + p1[i] * beta;
			}
			/*end*/
		}
		/* *********************************************
		* calculate norm of residual & check convergence
		* **********************************************/
		EXCHANGE(cF);
		res = getResidual(r * pC,cF,sync);
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
					if(gInterMesh[patchi].to == source) 
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
