#include <cuda.h>
#include "solve.h"

/*number of threads in a block*/
static const Int nThreads = 128;

/*Matrix vector multiply*/
template <class T>
__global__
void cudaMul(const Int* const rows,
             const Int* const cols,
             const Scalar* const an,
             const Int N,
             const T* const x, 
             T* y
             ) {
    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)  {
        const Int start = rows[i];
        const Int end = rows[i + 1];
        T res = an[start] * x[cols[start]];

        for (Int j = start + 1; j < end; j++)
            res -= an[j] * x[cols[j]];
        y[i] = res;
    }
}
/*jacobi solver*/
template<class T>
__global__
void cudaJacobi(const Int* const rows,
                 const Int* const cols,
                 const Scalar* const an,
                 const T* const cF,
                 T* const cF1,
                 const T* const Su,
                 T* r,
                 const Int N, 
                 Scalar omega
                 ) {
    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)  {
        const Int start = rows[i];
        const Int end = rows[i + 1];
        T res = Su[i], val = cF[i];

        for (Int j = start + 1; j < end; j++)
            res += an[j] * cF[cols[j]];
        res /= an[start];

        r[i] = -val;
        val *= (1 - omega);
        val += res * (omega);
        r[i] += val;
        cF1[i] = val;
    }
}
/*Taxpy*/
template<class T,class T1>
__global__
void cudaTaxpy(const Int N,
               const T1 alpha,
               const T* const x,
               const T* const y,
               T* const z
               ) {
    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)  {
        T temp;
        temp = x[i];
        temp *= alpha;
        temp += y[i];
        z[i] = temp;
    }
}
/*Txmy*/
template<class T,class T1>
__global__
void cudaTxmy(const Int N,
              const T* const x,
              const T1* const y,
              T* const z
              ) {
    Int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)  {
        T temp;
        temp = x[i];
        temp *= y[i];
        z[i] = temp;
    }
}
/*Tdot*/
template <class T>
__global__ 
void Tdot(const T* const a, 
          const T* const b, 
          T* const c, 
          const Int N
          ) {
    __shared__ T cache[nThreads];
    Int tid = threadIdx.x + blockIdx.x * blockDim.x;
    Int cacheIndex = threadIdx.x;

    T   temp = T(0),val;
    while (tid < N) {
        val = a[tid];
        val *= b[tid];
        temp += val;
        tid += blockDim.x * gridDim.x;
    }  
    cache[cacheIndex] = temp;
    
    __syncthreads();

    Int i = blockDim.x / 2;
    while (i != 0) {
        if (cacheIndex < i)
            cache[cacheIndex] += cache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        c[blockIdx.x] = cache[0];
}
template<class T>
__host__ 
T cudaTdot(T* x,
           T* y,
           T* d_sum,
           T* sum,
           const Int nBlocks32,
           const Int N
           ) {
    Tdot <<< nBlocks32, nThreads >>> (x,y,d_sum,N);
    cudaMemcpy(sum,d_sum,nBlocks32 * sizeof(T),cudaMemcpyDeviceToHost);
    T c = T(0);
    for (Int i = 0; i < nBlocks32; i++)
        c += sum[i];
    return c;
}
/* **************************************
 *   CSR - compressed sparse row format
 *       * Used for on GPU computation
 *       * Propably for AMG too
 * **************************************/
template <class T>
class CSRMatrix {
public:
    std::vector<Int>  rows;
    std::vector<Int>  cols;
    std::vector<Scalar> an;
    std::vector<Scalar> anT;
    std::vector<T> cF;
    std::vector<T> Su;
public:
    template <class T1>
    CSRMatrix(const MeshMatrix<T1>& A) {
        using namespace Mesh;
        using namespace DG;
        const Int N  = A.ap.size();
        const Int NN = A.ap.size() + 
                       A.an[0].size() + 
                       A.an[1].size() +
                       (NPMAT ? (A.adg.size() - A.ap.size()) : 0); 
        register Int i,f;

        /*resize*/
        cF.resize(N);
        Su.resize(N);
        rows.reserve(N + 1);
        cols.reserve(NN);
        an.reserve(NN);
        anT.reserve(NN);

        /*source term*/
        for(i = 0;i < N;i++) {
            Su[i] = A.Su[i];
            cF[i] = (*A.cF)[i];
        }

        /*fill matrix in CSR format.Diagonal element 
          is always at the start of a row */
        Int cn = 0;
        for(ii = 0;ii < gCells.size();ii++) {
            Cell& c = gCells[ii];
            for(Int j = 0;j < NP;j++) {
                Int i = ii * NP + j;

                rows.push_back(cn);

                an.push_back(A.ap[i]);
                anT.push_back(A.ap[i]);
                cols.push_back(i);
                cn++;

                forEach(c,k) {
                    f = c[k];
                    if(i == gFO[f]) {
                        an.push_back(A.an[1][f]);
                        anT.push_back(A.an[0][f]);
                        cols.push_back(gFN[f]);
                        cn++;
                    } else {
                        an.push_back(A.an[0][f]);
                        anT.push_back(A.an[1][f]);
                        cols.push_back(gFO[f]);
                        cn++;
                    }
                }
                
                if(NPMAT) {
                    for(Int k = 0;k < NP;k++) {
                        if(k == j) continue;
                        an.push_back(A.adg[ii * NPMAT + j * NP + k]);
                        anT.push_back(A.adg[ii * NPMAT + k * NP + j]);
                        cols.push_back();
                        cn++;
                    }
                }
            }
        }
        /*push extra row*/
        rows.push_back(cn);
    }
    /*IO*/
    friend std::ostream& operator << (std::ostream& os, const CSRMatrix& p) {
        os << p.rows << std::endl;
        os << p.cols << std::endl;
        os << p.an << std::endl;
        os << p.Su << std::endl;
        return os;
    }
    friend std::istream& operator >> (std::istream& is, CSRMatrix& p) {
        is >> p.rows;
        is >> p.cols;
        is >> p.an;
        is >> p.Su;
        return is;
    }
    /*end*/
};
/***********************************************
 * Template class to solve equations on GPU
 *      Solver must do many iterations to compensate
 *      for the latency caused by copying matrix
 *      from host to device.
 ***********************************************/
template<class T>
__host__
void SolveT(const MeshMatrix<T>& M) {
    const Int N = Mesh::gBCSfield;
    const Int Nall = M.ap.size();
    const Int nBlocks = (N + nThreads - 1) / nThreads;
    const Int nBlocks32 = ((nBlocks > 32) ? 32 : nBlocks);

    //info
    if(M.flags & M.SYMMETRIC)
        MP::printH("Symmetric  : ");
    else
        MP::printH("Asymmetric : ");
    if(Controls::Solver == Controls::SOR)
        MP::print("SOR :");
    else
        MP::print("PCG :");

    /*******************************
     *  variables on host & device
     *******************************/
    Int*   d_rows;
    Int*   d_cols;
    Scalar*  d_an;
    Scalar*  d_anT;
    Scalar*  d_pC;
    T*       d_cF;
    T*       d_Su;
    //PCG
    T*       d_r,*d_r1;
    T*       d_p,*d_p1,*d_AP,*d_AP1;
    T        alpha,beta,o_rr,oo_rr;
    T        local_res[2];
    //reduction
    T*       sum,*d_sum;

    /*********************************
     * allocate memory on device
     ********************************/
    {
        CSRMatrix<T> A(M);  
        cudaMalloc((void**) &d_rows,A.rows.size() * sizeof(Int));
        cudaMalloc((void**) &d_cols,A.cols.size() * sizeof(Int));
        cudaMalloc((void**) &d_an,  A.an.size() * sizeof(Scalar));
        cudaMalloc((void**) &d_cF,  Nall * sizeof(T));
        cudaMalloc((void**) &d_Su,  Nall * sizeof(T));

        cudaMemcpy(d_rows ,&A.rows[0] ,A.rows.size() * sizeof(Int),  cudaMemcpyHostToDevice);
        cudaMemcpy(d_cols ,&A.cols[0] ,A.cols.size() * sizeof(Int),  cudaMemcpyHostToDevice);
        cudaMemcpy(d_an   ,&A.an[0]   ,A.an.size() * sizeof(Scalar), cudaMemcpyHostToDevice);
        cudaMemcpy(d_cF   ,&A.cF[0]   ,Nall *   sizeof(T),    cudaMemcpyHostToDevice);
        cudaMemcpy(d_Su   ,&A.Su[0]   ,Nall *   sizeof(T),    cudaMemcpyHostToDevice);

        cudaMalloc((void**) &d_r, Nall * sizeof(T));
        cudaMalloc((void**) &d_sum, nBlocks32 * sizeof(T));
        sum = (T*) malloc(nBlocks32 * sizeof(T));

        if(Controls::Solver == Controls::SOR) {
            cudaMalloc((void**) &d_AP,Nall * sizeof(T));
            cudaMemcpy( d_AP,d_cF,Nall * sizeof(T),cudaMemcpyDeviceToDevice);
        } else if(Controls::Solver == Controls::PCG) {
            cudaMalloc((void**) &d_p,   Nall * sizeof(T));
            cudaMalloc((void**) &d_AP,  Nall * sizeof(T));
            {
                ScalarCellField pC = 1./M.ap;
                cudaMalloc((void**) &d_pC,N * sizeof(Scalar));
                cudaMemcpy(d_pC,&pC[0],N * sizeof(Scalar),cudaMemcpyHostToDevice);
            }
            if(!(M.flags & M.SYMMETRIC)) {
                cudaMalloc((void**) &d_r1,   Nall * sizeof(T));
                cudaMalloc((void**) &d_p1,   Nall * sizeof(T));
                cudaMalloc((void**) &d_AP1,  Nall * sizeof(T));
                cudaMalloc((void**) &d_anT,A.anT.size() * sizeof(Scalar));
                cudaMemcpy(d_anT,&A.anT[0],A.anT.size() * sizeof(Scalar), cudaMemcpyHostToDevice);
            }
        }
    }

    /*CG*/
    if(Controls::Solver == Controls::PCG) {
        cudaMemset(d_r,0,Nall * sizeof(T));
        cudaMemset(d_p,0,Nall * sizeof(T));
        cudaMul   <<< nBlocks, nThreads >>> (d_rows,d_cols,d_an,N,d_cF,d_AP);
        cudaTaxpy <<< nBlocks, nThreads >>> (N,Scalar(-1),d_AP,d_Su,d_r);
        cudaTxmy  <<< nBlocks, nThreads >>> (N,d_r,d_pC,d_p);
        o_rr = cudaTdot(d_r,d_p,d_sum,sum,nBlocks32,N);
    }
    /*BiCG*/
    if(!(M.flags & M.SYMMETRIC) && (Controls::Solver == Controls::PCG)) {
        cudaMemcpy(d_r1,d_r,Nall * sizeof(T), cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_p1,d_p,Nall * sizeof(T), cudaMemcpyDeviceToDevice);
    }
    //iterate until convergence
    Scalar res = 0;
    Int iterations = 0;

    /* **************************
     * Iterative solvers
     * *************************/
    while(iterations < Controls::max_iterations) {
        /*counter*/
        iterations++;

        /*select solver*/
        if(Controls::Solver == Controls::SOR) {
            iterations++;
            cudaJacobi <<< nBlocks, nThreads >>> (d_rows,d_cols,d_an,d_cF,d_AP,d_Su,d_r,N,Controls::SOR_omega);
            cudaJacobi <<< nBlocks, nThreads >>> (d_rows,d_cols,d_an,d_AP,d_cF,d_Su,d_r,N,Controls::SOR_omega);
        } else if(M.flags & M.SYMMETRIC) {
            /*conjugate gradient   : from wiki*/
            cudaMul <<< nBlocks, nThreads >>> (d_rows,d_cols,d_an,N,d_p,d_AP);
            oo_rr = cudaTdot(d_p,d_AP,d_sum,sum,nBlocks32,N);
            alpha = sdiv(o_rr , oo_rr);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,alpha,d_p,d_cF,d_cF);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,-alpha,d_AP,d_r,d_r);
            oo_rr = o_rr;
            cudaTxmy <<< nBlocks, nThreads >>> (N,d_r,d_pC,d_AP);
            o_rr = cudaTdot(d_r,d_AP,d_sum,sum,nBlocks32,N);
            beta = sdiv(o_rr , oo_rr);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,beta,d_p,d_AP,d_p);
            /*end*/
        } else {
            /* biconjugate gradient : from wiki */
            cudaMul <<< nBlocks, nThreads >>> (d_rows,d_cols,d_an,N,d_p,d_AP);
            cudaMul <<< nBlocks, nThreads >>> (d_rows,d_cols,d_anT,N,d_p1,d_AP1);
            oo_rr = cudaTdot(d_p1,d_AP,d_sum,sum,nBlocks32,N);
            alpha = sdiv(o_rr , oo_rr);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,alpha,d_p,d_cF,d_cF);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,-alpha,d_AP,d_r,d_r);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,-alpha,d_AP1,d_r1,d_r1);
            oo_rr = o_rr;
            cudaTxmy <<< nBlocks, nThreads >>> (N,d_r,d_pC,d_AP);
            cudaTxmy <<< nBlocks, nThreads >>> (N,d_r1,d_pC,d_AP1);
            o_rr = cudaTdot(d_r1,d_AP,d_sum,sum,nBlocks32,N);
            beta = sdiv(o_rr , oo_rr);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,beta,d_p,d_AP,d_p);
            cudaTaxpy <<< nBlocks, nThreads >>> (N,beta,d_p1,d_AP1,d_p1);
        }

        /* *********************************************
        * calculate norm of residual & check convergence
        * **********************************************/
        local_res[0] = cudaTdot(d_r,d_r,d_sum,sum,nBlocks32,N);
        local_res[1] = cudaTdot(d_cF,d_cF,d_sum,sum,nBlocks32,N);
        res = sqrt(mag(local_res[0]) / mag(local_res[1]));
        
        /*check convergence*/
        if(res <= Controls::tolerance)
            break;
    }

    /*****************************
     *  Copy result back to cpu
     *****************************/
    //copy result
    cudaMemcpy(&((*M.cF)[0]), d_cF, N * sizeof(T), cudaMemcpyDeviceToHost);

    //update boundary conditons
    updateExplicitBCs(*M.cF);

    //info
    MP::print("Iterations %d Residue: %.5e\n",iterations,res);
    /*********************************
     * free device memory
     ********************************/
    {
        cudaFree(d_rows);
        cudaFree(d_cols);
        cudaFree(d_an);
        cudaFree(d_cF);
        cudaFree(d_Su);

        cudaFree(d_r);
        cudaFree(d_sum);
        free(sum);

        if(Controls::Solver == Controls::SOR) {
            cudaFree(d_AP);
        } else if(Controls::Solver == Controls::PCG) {
            cudaFree(d_p);
            cudaFree(d_AP);
            cudaFree(d_pC);
            if(!(M.flags & M.SYMMETRIC)) {
                cudaFree(d_r1);
                cudaFree(d_p1);
                cudaFree(d_AP1);
                cudaFree(d_anT);
            }
        }
    }
    /******************
     *    END
     ******************/
}

/***************************
 * Explicit instantiations
 ***************************
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
