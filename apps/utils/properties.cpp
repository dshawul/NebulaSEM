#include "properties.h"

namespace Fluid {
    Scalar density;    /**< Density */
    Scalar viscosity;  /**< Viscosity */
    Scalar Pr;         /**< Prandtl number */
    Scalar Prt;        /**< Turbulent Prandtly number */
    Scalar T0;         /**< Reference temperature */
    Scalar P0;         /**< Reference pressure */
    Scalar cp;         /**< Specific heat at constant pressure */
    Scalar cv;         /**< Specific heat at constant volume */
    Scalar beta;       /**< Coefficient of thermal expansion (1/T0) */
    
    void enroll(Util::ParamList& params) {
        density = 1.177;
        viscosity = 1.568e-5;
        Pr = 0.9;
        Prt = 0.7;
        beta = 3.33e-3;
        T0 = 300;
        P0 = 101325;
        cp = 1004.67;
        cv = 715.5;
    
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
