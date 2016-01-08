#ifndef __MIXING_LENGTH_H
#define __MIXING_LENGTH_H

#include "turbulence.h"

struct MixingLength_Model : public EddyViscosity_Model {
    /*model coefficients*/
    Scalar mixingLength;
    Scalar C;
    Int wallDamping;

    /*mixing length field*/
    ScalarCellField lm;
    Scalar kappa;

    /*constructor*/
    MixingLength_Model(VectorCellField&,ScalarFacetField&,ScalarCellField&,Scalar&);

    /*others*/
    virtual void enroll();
    virtual void calcEddyViscosity(const TensorCellField& gradU);
    virtual void applyWallFunction(Int f,LawOfWall& low);
    virtual ScalarCellField getK();
    virtual void calcLengthScale() { 
        lm = mixingLength;
    }
};

#endif
