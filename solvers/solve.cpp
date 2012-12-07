#include "solve.h"

/* *********************************************************************
 *  Solve system of linear equations iteratively
 * *********************************************************************/
template<class type>
void SolveT(const MeshMatrix<type>& M) {

	using namespace Mesh;
	MeshField<type,CELL> r,p,AP;
	MeshField<type,CELL> r1(false),p1(false),AP1(false);   /* Allocate only if BiCG is used*/
	MeshField<type,CELL>& cF = *M.cF;
	MeshField<type,CELL>& buffer = AP;
	ScalarCellField pC = (1 / M.ap);           /* Jacobi preconditioner */
	Scalar res,ires;
	Int i,j,iterations = 0;
	type alpha,beta,o_rr = type(0),oo_rr;
	type local_res[2];

	/*for cluster use*/
	bool converged = false;
	bool init_exchange = true;
	bool print = (MP::host_id == 0);
	int  end_count = gInterMesh.size();
	std::vector<bool> sent_end(gInterMesh.size(),false);

	/*identify matrix solver type*/
	if(print) {
		if(M.flags & M.SYMMETRIC)
			MP::printH("Symmetric  : ");
		else
			MP::printH("Asymmetric : ");
		if(Controls::Solver == Controls::SOR)
			MP::print("SOR :");
		else
			MP::print("PCG :");
	}
	/*initial residual*/
	r = M.Su - M * cF;
	for(i = gBCellsStart;i < r.size();i++)
		r[i] = type(0);
	local_res[0] = type(0);
	local_res[1] = type(0); 
	for(i = 0;i < gBCellsStart;i++) {
		local_res[0] += (r[i] * r[i]);
		local_res[1] += (cF[i] * cF[i]);
	}
	res = ires = sqrt(mag(local_res[0]) / mag(local_res[1]));

	/*CG*/
	if(Controls::Solver == Controls::PCG) {
		for(i = gBCellsStart;i < r.size();i++)
			p[i] = type(0);
		o_rr = type(0);
		for(i = 0;i < gBCellsStart;i++) {
			p[i] = r[i] * pC[i];
			o_rr += r[i] * p[i];
		}
		/*BiCG*/
		if(!(M.flags & M.SYMMETRIC)) {
			r1.allocate();
			p1.allocate();
			AP1.allocate();
			r1 = r;
			p1 = p;
		}
	}
    /* **************************
	 * Iterative solvers
	 * *************************/
	while(iterations < Controls::max_iterations) {
		/*counter*/
		iterations++;

		/*select solver*/
		if(Controls::Solver == Controls::SOR) {
			/*default solver SOR*/
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

				r[i] = -cF[i];
				cF[i] = cF[i] * (1 - Controls::SOR_omega) + 
						ncF * (Controls::SOR_omega);
				r[i] += cF[i];
			}
			/*end*/
		} else if(M.flags & M.SYMMETRIC) {
			/*conjugate gradient   : from wiki*/
			AP = M * p;
			oo_rr = type(0);
			for(i = 0;i < gBCellsStart;i++)
				oo_rr += p[i] * AP[i];
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
			beta = sdiv(o_rr , oo_rr);
			for(i = 0;i < gBCellsStart;i++) {
				p[i] = AP[i] + p[i] * beta;
			}
			/*end*/
		} else {
			/* biconjugate gradient : from wiki */
			AP = M * p;
			AP1 = M ^ p1;
			oo_rr = type(0);
			for(i = 0;i < gBCellsStart;i++)
				oo_rr += p1[i] * AP[i];
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
		local_res[0] = type(0);
		local_res[1] = type(0); 
		for(i = 0;i < gBCellsStart;i++) {
			local_res[0] += (r[i] * r[i]);
			local_res[1] += (cF[i] * cF[i]);
		}
		res = sqrt(mag(local_res[0]) / mag(local_res[1]));
		
		if(res <= Controls::tolerance)
			converged = true;
PROBE:
		/* **********************************************************
		 * update inter mesh boundary conditons
		 *  -> Communication is NOT forced on every iteration,
		 *  rather a non-blocking probe is used to check for message.
		 ************************************************************/
		{
			int source,message_id,cn;

			/*initialize exchange*/
			if(init_exchange) {
				init_exchange = false;
				for(i = 0;i < gInterMesh.size();i++) {
					interBoundary& b = gInterMesh[i];
					if(b.from < b.to) {
						IntVector& f = *(b.f);
						for(j = 0;j < f.size();j++)
							buffer[j] = cF[gFN[f[j]]];
						MP::send(&buffer[0],f.size(),b.to,MP::FIELD);
					}
				}
			}
			/*probe*/
			while(MP::iprobe(source,message_id)) {
				/*find the boundary*/
				for(i = 0;i < gInterMesh.size();i++) {
					if(gInterMesh[i].to == source) break;
				}
				interBoundary& b = gInterMesh[i];
				/*parse message*/
				if(message_id == MP::FIELD) {
					IntVector& f = *(b.f);
					MP::recieve(&buffer[0],f.size(),source,message_id);
					type temp;

					/*residual*/
					Scalar res;
					local_res[0] = type(0);
					local_res[1] = type(0); 
					for(j = 0;j < f.size();j++) {
						cn = gFN[f[j]];
						temp = buffer[j] - cF[cn];
						local_res[0] += (temp * temp);
						local_res[1] += (cF[cn] * cF[cn]);
					}
					res = sqrt(mag(local_res[0]) / mag(local_res[1]));

					/*if change is small stop communication*/
					if(res > Controls::tolerance && iterations < Controls::max_iterations) {
						/*exchange and send back*/
						for(j = 0;j < f.size();j++) {
							cn = gFN[f[j]];
							temp = cF[cn];
							cF[cn] = buffer[j];
							buffer[j] = temp;
						}
						MP::send(&buffer[0],f.size(),source,message_id);
						/*set flag off*/
						converged = false;
					} else {
						if(!sent_end[i]) {
							MP::send(source,MP::END);
							sent_end[i] = true;
						}
					}
				} else if(message_id == MP::END) {
					MP::recieve(source,message_id);
					end_count--;
					if(!sent_end[i]) {
						MP::send(source,MP::END);
						sent_end[i] = true;
					}
				}
			}
		}
		/* **************************************
		* Continue iteration if we have updates
		* **************************************/
		if(iterations == Controls::max_iterations)
			converged = true;
		if(converged) {
			if(end_count > 0) goto PROBE;
			else break;
		}
	}

	/*solver info*/
	if(print)
		MP::print("Iterations %d Initial Residual %.5e Final Residual %.5e\n",iterations,ires,res);

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
