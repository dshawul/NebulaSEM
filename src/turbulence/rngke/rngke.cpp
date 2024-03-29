#include "rngke.h"
/**
RNG-KE turbulence model constructor
References:
    http://www.cfd-online.com/Wiki/RNG_k-epsilon_model
 */
RNG_KE_Model::RNG_KE_Model(VectorCellField& tU,VectorCellField& tFc,ScalarFacetField& tF,ScalarCellField& trho,ScalarCellField& tmu) :
    KE_Model(tU,tFc,tF,trho,tmu),
    eta0(4.38),
    beta(0.012)
{
    Cmu = 0.0845;
    SigmaK = 0.7194;
    SigmaX = 0.7194;
    C1x = 1.42;
    C2x = 1.68;
}
/** Enroll parameters */
void RNG_KE_Model::enroll() {
    using namespace Util;
    KE_Model::enroll();
    params.enroll("eta0",&eta0);
    params.enroll("beta",&beta);
}
/** Calculate eddy viscosity */
void RNG_KE_Model::calcEddyViscosity(const TensorCellField& gradU) {
    /*calculate C2eStar*/
    {
        ScalarCellField eta = sqrt(getS2(gradU)) * (k / x);
        C2eStar = C2x + Cmu * pow(eta,3.0) * (1 - eta / eta0) / 
            (1 + beta * pow(eta,3.0));
        C2eStar = max(C2eStar,Scalar(0.0));
    }
    /*calculate viscosity*/
    KE_Model::calcEddyViscosity(gradU);
}
/** Solver transport of turbulent dissipation and kinetic energy */
void RNG_KE_Model::solve() {
    ScalarCellMatrix M;
    ScalarCellField eff_mu;

    /*turbulent dissipation*/
    eff_mu = eddy_mu / SigmaX + mu;
    M = transport<Scalar>(x, Fc, F, eff_mu, x_UR,
            (C1x * Pk * x / k),
            -(C2eStar * rho * x / k), 0, &rho);
    FixNearWallValues(M);
    Solve(M);
    x = max(x,Constants::MachineEpsilon);

    /*turbulent kinetic energy*/
    eff_mu = eddy_mu / SigmaK + mu;
    M = transport<Scalar>(k, Fc, F, eff_mu, k_UR,
            Pk,
            -(rho * x / k), 0, &rho);
    if(wallModel == STANDARD)
        FixNearWallValues(M);
    Solve(M);
    k = max(k,Constants::MachineEpsilon);
}
