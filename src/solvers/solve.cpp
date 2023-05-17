#include "solve.h"

/* *********************************************************************
 *  Solve system of linear equations iteratively
 * *********************************************************************/

/**
  Calculate global residual
 */
template<class type, ENTITY entity>
Scalar getResidual(const MeshField<type,entity>& r,
        const MeshField<type,entity>& cF,
        bool sync) {
    type res0, res1;
    res0 = type(0);
    res1 = type(0); 
    #pragma omp parallel for reduction(+:res0,res1)
    for(Int i = 0;i < Mesh::gBCSfield;i++) {
        res0 += (r[i] * r[i]);
        res1 += (cF[i] * cF[i]);
    }
    if(sync) {
        type global_res0, global_res1;
        MP::allreduce(&res0,&global_res0,1,MP::OP_SUM);
        MP::allreduce(&res1,&global_res1,1,MP::OP_SUM);
        res0 = global_res0;
        res1 = global_res1;
    }
    return sqrt(sdiv(mag(res0), mag(res1)));
}
/**
  Solve a system of linear equations Ax=B
 */
template<class T1, class T2, class T3>
void SolveT(const MeshMatrix<T1,T2,T3>& M) {
    using namespace Mesh;
    using namespace DG;
    MeshField<T3,CELL> r,p,AP = T3(0);
    MeshField<T3,CELL> r1(false),p1(false),AP1(false);   
    MeshField<T1,CELL>& cF = *M.cF;
    MeshField<T3,CELL>& buffer = AP;
    MeshField<T2,CELL> D = M.ap,iD = (T2(1) / M.ap);
    Scalar res,ires;
    T1 alpha,beta,o_rr = T1(0),oo_rr;
    Int iterations = 0;
    bool converged = false;

    /****************************
     * Parallel controls
     ***************************/
    int  end_count = 0;
    bool sync = (Controls::parallel_method == Controls::BLOCKED)
        && gInterMesh.size();
    std::vector<bool> sent_end(gInterMesh.size(),false);

    /****************************
     * Jacobi sweep
     ***************************/
#define JacobiSweep() {                             \
    cF = iD * getRHS(M,sync);                       \
}
    /****************************
     *  Forward/backward GS sweeps
     ****************************/
#define Sweep_(X,B,ci) {                            \
    Cell& c = gCells[ci];                           \
    forEachLgl(ii,jj,kk) {                          \
        Int index1 = INDEX4(ci,ii,jj,kk);           \
        T3 ncF = B[index1];                         \
        if(NPMAT) {                                 \
            T3 val(Scalar(0));                      \
            forEachLglX(i) {                                                \
                Int index2 = INDEX4(ci,i,jj,kk);                            \
                Int indexm = ci * NPMAT + INDEX_X(ii,jj,kk,i);              \
                val += X[index2] * M.adg[indexm];                           \
            }                                                               \
            forEachLglY(j) if(j != jj) {                                    \
                Int index2 = INDEX4(ci,ii,j,kk);                            \
                Int indexm = ci * NPMAT + INDEX_Y(ii,jj,kk,j);              \
                val += X[index2] * M.adg[indexm];                           \
            }                                                               \
            forEachLglZ(k) if(k != kk) {                                    \
                Int index2 = INDEX4(ci,ii,jj,k);                            \
                Int indexm = ci * NPMAT + INDEX_Z(ii,jj,kk,k);              \
                val += X[index2] * M.adg[indexm];                           \
            }                                                               \
            ncF += val;                             \
        }                                           \
        forEach(c,j) {                              \
            Int faceid = c[j];                      \
            for(Int n = 0; n < NPF;n++) {           \
                Int k = faceid * NPF + n;           \
                Int c1 = FO[k];                     \
                Int c2 = FN[k];                     \
                if(index1 == c1)                    \
                    ncF += X[c2] * M.ann[k];        \
                else if(index1 == c2)               \
                    ncF += X[c1] * M.ano[k];        \
            }                                       \
        }                                           \
        ncF *= iD[index1];                                  \
        X[index1] = X[index1] * (1 - Controls::SOR_omega) + \
            ncF * (Controls::SOR_omega);                    \
    }                                                       \
}
#define ForwardSweep(X,B) {                         \
    ASYNC_COMM<T1> comm(&X[0]);                     \
    comm.send();                                    \
    _Pragma("omp parallel for")                     \
    for(Int ci = 0;ci < gBCSI;ci++)                 \
        Sweep_(X,B,ci);                             \
    comm.recv();                                    \
    _Pragma("omp parallel for")                     \
    for(Int ci = gBCSI;ci < gBCS;ci++)              \
        Sweep_(X,B,ci);                             \
}
    /***********************************
     *  Forward/backward substitution
     ***********************************/
#define Substitute_(X,B,ci,forw,tr) {           \
    const Int index1 = INDEX4(ci,ii,jj,kk);     \
    T3 ncF = B[index1];                         \
    if(NPMAT) {                                 \
        T3 val(Scalar(0));                      \
        forEachLglX(i) {                                                        \
            Int index2 = INDEX4(ci,i,jj,kk);                                    \
            if((forw && (index2 < index1)) ||   (!forw && (index1 < index2))) { \
                Int indexm = ci * NPMAT +                                       \
                    (tr ? INDEX_TX(ii,jj,kk,i) : INDEX_X(ii,jj,kk,i));          \
                val += X[index2] * M.adg[indexm];                               \
            }                                                                   \
        }                                                                       \
        forEachLglY(j) if(j != jj) {                                            \
            Int index2 = INDEX4(ci,ii,j,kk);                                    \
            if((forw && (index2 < index1)) ||   (!forw && (index1 < index2))) { \
                Int indexm = ci * NPMAT +                                       \
                (tr ? INDEX_TY(ii,jj,kk,j) : INDEX_Y(ii,jj,kk,j));              \
                val += X[index2] * M.adg[indexm];                               \
            }                                                                   \
        }                                                                       \
        forEachLglZ(k) if(k != kk) {                                            \
            Int index2 = INDEX4(ci,ii,jj,k);                                    \
            if((forw && (index2 < index1)) ||   (!forw && (index1 < index2))) { \
                Int indexm = ci * NPMAT +                                       \
                (tr ? INDEX_TZ(ii,jj,kk,k) : INDEX_Z(ii,jj,kk,k));              \
                val += X[index2] * M.adg[indexm];                               \
            }                                                                   \
        }                                                                       \
        ncF += val;                                 \
    }                                               \
    if(isBoundary(ii,jj,kk)) {                      \
        forEach(c,j) {                              \
            Int faceid = c[j];                      \
            for(Int n = 0; n < NPF;n++) {           \
                Int k = faceid * NPF + n;           \
                Int c1 = FO[k];                     \
                Int c2 = FN[k];                     \
                if(index1 == c1) {                  \
                    if((forw && (c2 < c1)) ||       \
                      (!forw && (c1 < c2))) {       \
                        if(tr == 0)                 \
                           ncF += X[c2] * M.ann[k]; \
                        else                        \
                           ncF += X[c2] * M.ano[k]; \
                    }                               \
                } else if(index1 == c2) {           \
                    if((forw && (c2 > c1)) ||       \
                      (!forw && (c1 > c2)))         \
                        if(tr == 0)                 \
                          ncF += X[c1] * M.ano[k];  \
                        else                        \
                          ncF += X[c1] * M.ann[k];  \
                }                                   \
            }                                       \
        }                                           \
    }                                               \
    ncF *= iD[index1];                              \
    X[index1] = ncF;                                \
}
#define ForwardSub(X,B,TR) {                        \
    _Pragma("omp parallel for")                     \
    for(Int ci = 0;ci < gBCS;ci++)  {               \
        Cell& c = gCells[ci];                       \
        forEachLgl(ii,jj,kk)                        \
            Substitute_(X,B,ci,true,TR);            \
    }                                               \
}
#define BackwardSub(X,B,TR) {                       \
    _Pragma("omp parallel for")                     \
    for(int ci = gBCS - 1; ci >= 0; ci--)    {      \
        Cell& c = gCells[ci];                       \
        forEachLglR(ii,jj,kk)                       \
            Substitute_(X,B,ci,false,TR);           \
    }                                               \
}
#define DiagSub(X,B) {                              \
    _Pragma("omp parallel for")                     \
    for(Int i = 0;i < gBCSfield;i++)                \
        X[i] = B[i] * iD[i];                        \
}
    /***********************************
     *  Preconditioners
     ***********************************/
#define precondition_(R,Z,TR) {                     \
    using namespace Controls;                       \
    if(Preconditioner == Controls::NOPR) {          \
        Z = R;                                      \
    } else if(Preconditioner == Controls::DIAG) {   \
        DiagSub(Z,R);                               \
    } else {                                        \
        if(Controls::Solver == Controls::PCG) {     \
            ForwardSub(Z,R,TR);                     \
            Z = Z * D;                              \
            BackwardSub(Z,Z,TR);                    \
        }                                           \
    }                                               \
}
#define precondition(R,Z) precondition_(R,Z,0)
#define preconditionT(R,Z) precondition_(R,Z,1)
    /***********************************
     *  SAXPY and DOT operations
     ***********************************/
#define Taxpy(Y,I,X,alpha_) {                       \
    _Pragma("omp parallel for")                     \
    for(Int i = 0;i < gBCSfield;i++)                \
        Y[i] = I[i] + X[i] * alpha_;                \
}
#define Tdot(X,Y,sumt) {                            \
    T3 sum = T3(0);                                 \
    _Pragma("omp parallel for reduction(+:sum)")    \
    for(Int i = 0;i < gBCSfield;i++)                \
        sum += X[i] * Y[i];                         \
    sumt = sum;                                     \
}
    /***********************************
     *  Synchronized sum
     ***********************************/
#define REDUCE(typ,var) if(sync) {                  \
    typ t;                                          \
    MP::allreduce(&var,&t,1,MP::OP_SUM);            \
    var = t;                                        \
}
    /***********************************
     *  Residual
     ***********************************/
#define CALC_RESID() {                              \
    r = M.Su - mul(M,cF);                           \
    _Pragma("omp parallel for")                     \
    forEachS(r,k,gBCSfield)                         \
        r[k] = T3(0);                               \
    precondition(r,AP);                             \
    _Pragma("omp parallel for")                     \
    forEachS(AP,k,gBCSfield)                        \
        AP[k] = T3(0);                              \
    res = getResidual(AP,cF,sync);                  \
    if(Controls::Solver == Controls::PCG) {         \
        Tdot(r,AP,o_rr);                            \
        REDUCE(T1,o_rr);                            \
        p = AP;                                     \
        if(!(M.flags & M.SYMMETRIC)) {              \
            r1 = r;                                 \
            p1 = p;                                 \
        }                                           \
    }                                               \
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
            if(Controls::Preconditioner == Controls::SSOR) {
                /*SSOR pre-conditioner*/
                iD *= Controls::SOR_omega;
                D *=  (2.0 / Controls::SOR_omega - 1.0);    
            } else if(Controls::Preconditioner == Controls::DILU) {
                /*D-ILU(0) pre-conditioner*/
                #pragma omp parallel for
                for(Int ci = 0;ci < gBCS;ci++) {
                    Cell& c = gCells[ci];
                    forEachLgl(ii,jj,kk) {
                        Int index1 = INDEX4(ci,ii,jj,kk);
                        if(NPMAT) {
                            T2 val = T2(0);
                            forEachLglX(i) {
                                Int index2 = INDEX4(ci,i,jj,kk);
                                if(index1 > index2) {
                                    val += iD[index2] * 
                                        M.adg[ci * NPMAT + INDEX_X(ii,jj,kk,i)] *
                                        M.adg[ci * NPMAT + INDEX_TX(ii,jj,kk,i)];
                                }
                            }
                            forEachLglY(j) if(j != jj) {
                                Int index2 = INDEX4(ci,ii,j,kk);
                                if(index1 > index2) {
                                    val += iD[index2] * 
                                        M.adg[ci * NPMAT + INDEX_Y(ii,jj,kk,j)] *
                                        M.adg[ci * NPMAT + INDEX_TY(ii,jj,kk,j)];
                                }
                            }
                            forEachLglZ(k) if(k != kk) {
                                Int index2 = INDEX4(ci,ii,jj,k);
                                if(index1 > index2) {
                                    val += iD[index2] * 
                                        M.adg[ci * NPMAT + INDEX_Z(ii,jj,kk,k)] *
                                        M.adg[ci * NPMAT + INDEX_TZ(ii,jj,kk,k)];
                                }
                            }
                            D[index1] -= val;
                        }   
                        if(isBoundary(ii,jj,kk)) {
                            forEach(c,j) {                              
                                Int faceid = c[j];
                                for(Int n = 0; n < NPF;n++) {
                                    Int k = faceid * NPF + n;                           
                                    Int c1 = FO[k];                     
                                    Int c2 = FN[k];                     
                                    if(index1 == c1) {
                                        if(c2 > c1) D[c2] -= 
                                            (M.ano[k] * M.ann[k] * iD[c1]); 
                                    } else if(index1 == c2) {
                                        if(c1 > c2) D[c1] -= 
                                            (M.ano[k] * M.ann[k] * iD[c2]);     
                                    }       
                                }                               
                            }
                        }
                    }           
                }
                iD = (T2(1) / D);
            }
            /*end*/
        }
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
                Int buf_size = f.size() * NPF;
                /*send*/
                forEach(f,j) {
                    Int faceid = f[j];
                    Int offset = j * NPF;
                    for(Int n = 0; n < NPF;n++) {
                        Int k = faceid * NPF + n;
                        buffer[offset + n] = cF[FO[k]];
                    }
                }
                MP::send(&buffer[0],buf_size,b.to,MP::FIELD);
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
        if(Controls::Solver != Controls::PCG) {
            p = cF;
            /*Jacobi and SOR solvers*/
            if(Controls::Solver == Controls::JACOBI) {
                JacobiSweep();
            } else {
                ForwardSweep(cF,M.Su);
            }
            /*residual*/
            #pragma omp parallel for
            for(Int i = 0;i < gBCSfield;i++)
                AP[i] = cF[i] - p[i];
        } else if(M.flags & M.SYMMETRIC) {
            /*conjugate gradient*/
            AP = mul(M,p,sync);
            Tdot(p,AP,oo_rr);
            REDUCE(T1,oo_rr);
            alpha = sdiv(o_rr , oo_rr);
            Taxpy(cF,cF,p,alpha);
            Taxpy(r,r,AP,-alpha);
            precondition(r,AP);
            oo_rr = o_rr;
            Tdot(r,AP,o_rr);
            REDUCE(T1,o_rr);
            beta = sdiv(o_rr , oo_rr);
            Taxpy(p,AP,p,beta);
            /*end*/
        } else {
            /* biconjugate gradient*/
            AP = mul(M,p,sync);
            AP1 = mult(M,p1,sync);
            Tdot(p1,AP,oo_rr);
            REDUCE(T1,oo_rr);
            alpha = sdiv(o_rr , oo_rr);
            Taxpy(cF,cF,p,alpha);
            Taxpy(r,r,AP,-alpha);
            Taxpy(r1,r1,AP1,-alpha);
            precondition(r,AP);
            preconditionT(r1,AP1);
            oo_rr = o_rr;
            Tdot(r1,AP,o_rr);
            REDUCE(T1,o_rr);
            beta = sdiv(o_rr , oo_rr);
            Taxpy(p,AP,p,beta);
            Taxpy(p1,AP1,p1,beta);
            /*end*/
        }
        /* *********************************************
         * calculate norm of residual & check convergence
         * **********************************************/
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
        if(!sync) {
            int source,message_id;
            /*probe*/
            while(MP::iprobe(source,message_id,MP::FIELD)
                    || MP::iprobe(source,message_id,MP::END)) {
                /*find the boundary*/
                Int patchi;
                for(patchi = 0;patchi < gInterMesh.size();patchi++) {
                    if(gInterMesh[patchi].to == (Int)source) 
                        break;
                }
                /*parse message*/
                if(message_id == MP::FIELD) {
                    interBoundary& b = gInterMesh[patchi];
                    IntVector& f = *(b.f);
                    Int buf_size = f.size() * NPF;
    
                    /*recieve*/
                    MP::recieve(&buffer[0],buf_size,source,message_id);
                    forEach(f,j) {
                        Int faceid = f[j];
                        Int offset = j * NPF;
                        for(Int n = 0; n < NPF;n++) {
                            Int k = faceid * NPF + n;
                            cF[FN[k]] = buffer[offset + n];
                        }
                    }
    
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
                        continue;
                    }
    
                    /*send*/
                    forEach(f,j) {
                        Int faceid = f[j];
                        Int offset = j * NPF;
                        for(Int n = 0; n < NPF;n++) {
                            Int k = faceid * NPF + n;
                            buffer[offset + n] = cF[FO[k]];
                        }
                    }
                    MP::send(&buffer[0],buf_size,source,message_id);
    
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
    /****************************
     * Iteration info
     ***************************/
    if(MP::printOn) {
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
                case Controls::NOPR: MP::print("NONE-PCG :"); break;
                case Controls::DIAG: MP::print("DIAG-PCG :"); break;
                case Controls::SSOR: MP::print("SSOR-PCG :"); break;
                case Controls::DILU: MP::print("DILU-PCG :"); break;
            }
        }
        MP::print("Iterations %d Initial Residual "
                "%.5e Final Residual %.5e\n",iterations,ires,res);
    }
}
/**
  Solve a diagonal system
 */
template<class T1,class T2,class T3>
void SolveTexplicit(const MeshMatrix<T1,T2,T3>& M) {
    *M.cF = M.Su / M.ap;
    if(MP::printOn) {
        MP::printH("DIAG-DIAG:");
        MP::print("Iterations %d Initial Residual "
                "%.5e Final Residual %.5e\n",1,0.0,0.0);
    }
}
/***************************
 * Explicit instantiations
 ***************************/
#define SOLVE() {                           \
    applyImplicitBCs(A);                    \
    if(A.flags & A.DIAGONAL)                \
        SolveTexplicit(A);                  \
    else                                    \
        SolveT(A);                          \
    applyExplicitBCs(*A.cF,true,false);     \
}
void Solve(const MeshMatrix<Scalar>& A) {
    SOLVE();
}
void Solve(const MeshMatrix<Vector>& A) {
    SOLVE();
}
void Solve(const MeshMatrix<STensor>& A) {
    SOLVE();
}
void Solve(const MeshMatrix<Tensor>& A) {
    SOLVE();
}
#undef SOLVE
/* ********************
 *        End
 * ********************/
