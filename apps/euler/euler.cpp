#include "solve.h"
#include "iteration.h"
#include "wrapper.h"
#include "properties.h"


/**
  Explicit gradient operator
 */
#define GRADD(im,jm,km) {                           \
    Int index1 = INDEX4(ci,im,jm,km);               \
    Vector dpsi_ij;                                 \
    DPSI(dpsi_ij,im,jm,km);                         \
    dpsi_ij = dot(Jin,dpsi_ij);                     \
    r[index1] -= mul(dpsi_ij,p[index]);             \
}

#define GRAD(T1,T2)                                                                             \
    inline auto mgradf(const MeshField<T2,CELL>& p) {                                            \
        using namespace Mesh;                                                                   \
        using namespace DG;                                                                     \
        MeshField<T1,CELL> r;                                                                   \
                                                                                                \
        r = sum(mul(fN,cds(p)));                                                                \
                                                                                                \
        /*if(NPMAT) {                                                                             \
            _Pragma("omp parallel for")                                                         \
            _Pragma("acc parallel loop copyin(p,gBCS)")                                         \
            for(Int ci = 0; ci < gBCS;ci++) {                                                   \
                forEachLgl(ii,jj,kk) {                                                          \
                    Int index = INDEX4(ci,ii,jj,kk);                                            \
                    Tensor Jin = Jinv[index] * cV[index];                                       \
                    forEachLglX(i) GRADD(i,jj,kk);                                              \
                    forEachLglY(j) if(j != jj) GRADD(ii,j,kk);                                  \
                    forEachLglZ(k) if(k != kk) GRADD(ii,jj,k);                                  \
                }                                                                               \
            }                                                                                   \
        }                                                                                       \
                                                                                                \
        fillBCs(r,false,p.fIndex);*/                                                              \
                                                                                                \
        return r;                                                                               \
    }

GRAD(Vector,Scalar);
GRAD(Tensor,Vector);
#undef GRAD
#undef GRADD

#define mgradi(x) (mgradf(x)  / Mesh::cV)


/**
  \verbatim
  Euler equations solver
  ~~~~~~~~~~~~~~~~~~~~~~
  d(rho*(1))/dt + div(1,F,0) = 0
  d(rho*(U))/dt + div(U,F,0) = -grad(p) + rho * gravity
  d(rho*(T))/dt + div(T,F,0) = 0
  \endverbatim
 */
void euler(std::istream& input) {
    /*Solver specific parameters*/
    Scalar pressure_UR = Scalar(0.5);
    Scalar velocity_UR = Scalar(0.8);
    Scalar t_UR = Scalar(0.8);
    Int buoyancy = 1;
    Int diffusion = 1;

    /*fluid properties*/
    {
        Util::ParamList params("general");
        Fluid::enroll(params);
        params.read(input);
    }

    /*euler options*/
    Util::ParamList params("euler");
    params.enroll("velocity_UR", &velocity_UR);
    params.enroll("pressure_UR", &pressure_UR);
    params.enroll("t_UR", &t_UR);
    Util::Option* op;
    op = new Util::BoolOption(&buoyancy);
    params.enroll("buoyancy",op);
    op = new Util::BoolOption(&diffusion);
    params.enroll("diffusion",op);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {

        ScalarCellField p("p",READWRITE);
        VectorCellField U("U",READWRITE);
        ScalarCellField T("T",READWRITE);
        ScalarCellField rho("rho",WRITE);
        VectorCellField g(false);

        /*Read fields*/
        Iteration it(ait.get_step());
        VectorCellField Fc;
        ScalarFacetField F;
        ScalarCellField rhof;
        ScalarCellField mu;

        /*gas constants*/
        using namespace Fluid;
        Scalar p_gamma, R, psi, iPr;
        R = cp - cv;
        p_gamma = cp / cv;
        psi = 1 / (R * T0);
        iPr = 1 / Pr;

        /*buoyancy*/
        ScalarCellField p_ref, rho_ref;
        if(buoyancy) {
            ScalarCellField gh;

            /*gravity vector*/
            g.construct("gravity",WRITE);
            if(Controls::is_spherical) {
                g = -unit(Mesh::cC) * mag(Controls::gravity);
                gh = -(mag(Mesh::cC) - Controls::sphere_radius) * mag(Controls::gravity);
            } else {
                g = Controls::gravity;
                gh = dot(g,Mesh::cC);
            }
            Mesh::fixedBCs<Vector>(U,g);
            applyExplicitBCs(g,true);
            g.write(0);

            /*compute hydrostatic pressure*/
            p_ref = P0 * pow(1.0 + gh / (cp * T), cp / R);
            Mesh::scaleBCs<Scalar>(p,p_ref,1.0);
            applyExplicitBCs(p_ref,true);

        } else {
            p_ref = P0;
        }
        p += p_ref;

        /*calculate rho*/
        rho_ref = (P0 / (R*T0)) * pow(p_ref / P0, 1 / p_gamma);
        Mesh::scaleBCs<Scalar>(p,rho_ref,psi);
        applyExplicitBCs(rho_ref,true);

        rho = (P0 / (R*T)) * pow(p / P0, 1 / p_gamma);
        Mesh::scaleBCs<Scalar>(p,rho,psi);
        applyExplicitBCs(rho,true);

        rho.write(0);

        /*Time loop*/
        for (; !it.end(); it.next()) {
            /*fluxes*/
            rhof = rho;
            Fc = flxc(rho * U);
            F = flx(rho * U);

            /*artificial viscosity*/
            if(diffusion) mu = rho * viscosity;
            else mu = Scalar(0);

            /*rho-equation*
            {
                ScalarCellMatrix M;
                M = convection(rho, flxc(U), flx(U), pressure_UR);
                Solve(M);
            }
            {
                ScalarCellField rhod = fabs(rho - rho_ref);
                Scalar maxv = reduce_max(rhod);
                Scalar minv = reduce_min(rhod);
                Scalar avgv = reduce_sum(rhod) / rhod.size();
                MP::printH("rho: Max: %g Min: %g Avg: %g\n",maxv,minv,avgv);
            }
            /*T-equation*/
            using namespace Util;
            ScalarFacetField Tf = cds(T);
            VectorCellField sN = sum(Mesh::fN) / Mesh::cV;
            std::cout << sN << std::endl;
            exit(0);
            VectorCellField gT = mgradi(T);
            for(Int i  = 0; i < Mesh::gBCSfield; i++)
                std::cout << i << ". " << T[i] << " " << gT[i] << std::endl;
            {
                ScalarCellMatrix M;
                //M = transport(T, Fc, F, mu * iPr, t_UR, &rho, &rhof);
                ScalarCellField mmu = mu * iPr;
                //M = div(T,Fc,F,&mmu) - lapf(T,mmu);
                M = -lap(T,mmu,true);
                addTemporal<1>(M,t_UR,&rho,&rhof);
                Solve(M);
            }

            /*calculate p*/
            p = P0 * pow((rho*T*R) / P0, p_gamma);
            applyExplicitBCs(p,true);
            {
                ScalarCellField rhod = fabs(p - p_ref);
                Scalar maxv = reduce_max(rhod);
                Scalar minv = reduce_min(rhod);
                Scalar avgv = reduce_sum(rhod) / rhod.size();
                MP::printH("pd: Max: %g Min: %g Avg: %g\n",maxv,minv,avgv);
            }

            /*U-equation*/
            {
                VectorCellField Sc = -gradi(p - p_ref);
                /*buoyancy*/
                if(buoyancy)
                    Sc += (rho - rho_ref) * g;
                {
                    Scalar maxv = reduce_max(Sc);
                    Scalar minv = reduce_min(Sc);
                    Vector avgv = reduce_sum(Sc) / Sc.size();
                    MP::printH("U src: Max: %g Min: %g Avg: %g\n",maxv,minv,mag(avgv));
                }
                /*solve*/
                VectorCellMatrix M;
                M = transport(U, Fc, F, mu, velocity_UR, Sc, Scalar(0), &rho, &rhof);
                Solve(M);
            }

            /*courant number*/
            if(Controls::state != Controls::STEADY)
                Mesh::calc_courant(U,Controls::dt);
            break;
        }
    }
}

/**
  \verbatim
  Main application entry point for euler solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   euler(Solver::input);
   Solver::Finalize();
   return 0;
}
