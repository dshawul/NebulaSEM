#ifndef __PROPERTIES_
#define __PROPERTIES_

#include "util.h"

/**
  general fluid properties
 */
namespace Fluid {

    extern Scalar density;
    extern Scalar viscosity;
    extern Scalar Pr;
    extern Scalar Prt;
    extern Scalar beta;
    extern Scalar T0;
    extern Scalar P0;
    extern Scalar cp;
    extern Scalar cv;

    void enroll(Util::ParamList& params);
}

#endif
