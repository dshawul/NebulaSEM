#include "properties.h"

namespace Fluid {
    Scalar density;
    Scalar viscosity;
    Scalar Pr;
    Scalar Prt;
    Scalar beta;
    Scalar T0;
    Scalar P0;
    Scalar cp;
    Scalar cv;
    
    void enroll(Util::ParamList& params) {
        density = 1.177;
        viscosity = 1.568e-5;
        Pr = 0.9;
        Prt = 0.7;
        beta = 3.33e-3;
        T0 = 300;
        P0 = 101325;
        cp = 1004.67;
    
        params.enroll("rho", &density);
        params.enroll("viscosity", &viscosity);
        params.enroll("Pr", &Pr);
        params.enroll("Prt", &Prt);
        params.enroll("beta", &beta);
        params.enroll("T0", &T0);
        params.enroll("P0", &P0);
        params.enroll("cp", &cp);
        params.enroll("cv", &cv);
    }
}
