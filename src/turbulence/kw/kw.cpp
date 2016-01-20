#include "kw.h"
/*
References:
    http://www.cfd-online.com/Wiki/Wilcox%27s_k-omega_model
*/
KW_Model::KW_Model(VectorCellField& tU,ScalarFacetField& tF,ScalarCellField& trho,ScalarCellField& tmu) :
    KX_Model(tU,tF,trho,tmu,"w")
{
    Cmu = 0.09;
    SigmaK = 2;
    SigmaX = 2;
    C1x = 5./9;
    C2x = 3./40;
}
void KW_Model::enroll() {
    using namespace Util;
    KX_Model::enroll();
    params.enroll("Cmu",&Cmu);
    params.enroll("SigmaK",&SigmaK);
    params.enroll("SigmaW",&SigmaX);
    params.enroll("C1w",&C1x);
    params.enroll("C2w",&C2x);
}
void KW_Model::solve() {
    ScalarCellMatrix M;
    ScalarCellField eff_mu;

    /*turbulent dissipation*/
    eff_mu = eddy_mu / SigmaX + mu;
    M = transport(x, U, F, eff_mu, x_UR,
                (C1x * Pk * x / k),
                -(C2x * rho * x), &rho);
    FixNearWallValues(M);
    Solve(M);
    x = max(x,Constants::MachineEpsilon);

    /*turbulent kinetic energy*/
    eff_mu = eddy_mu / SigmaK + mu;
    M = transport(k, U, F, eff_mu, k_UR,
                    Pk,
                    -(Cmu * rho * x), &rho);
    if(wallModel == STANDARD)
        FixNearWallValues(M);
    Solve(M);
    k = max(k,Constants::MachineEpsilon);
}
