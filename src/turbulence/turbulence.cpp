#include "turbulence.h"
#include "ke.h"
#include "rngke.h"
#include "realizableke.h"
#include "kw.h"
#include "les.h"

Int Turbulence_Model::turb_model = 2;
bool Turbulence_Model::bneedWallDist = false;

/** Register type of turbulence model desired */
void Turbulence_Model::RegisterTable(Util::ParamList& params) {
    Util::Option* op;
    op = new Util::Option(&turb_model,7,
        "NONE","MIXING_LENGTH","KE","RNG_KE","REALIZABLE_KE","KW","LES");
    params.enroll("turbulence_model",op);
}

/** Select type of turbulence model at run-time */
Turbulence_Model* Turbulence_Model::Select(VectorCellField& U,ScalarFacetField& F,
                    ScalarCellField& rho,ScalarCellField& mu) {
    /*turbulence model*/
    enum TurbModel {
        NONE,MIXING_LENGTH,KE,RNG_KE,REALIZABLE_KE,KW,LES
    };
    bneedWallDist = false;
    
    /*Select turbulence model*/
    Turbulence_Model* turb;
    switch(turb_model) {
        case KE:   
            turb = new KE_Model(U,F,rho,mu); 
            break;
        case MIXING_LENGTH:   
            bneedWallDist = true;
            turb = new MixingLength_Model(U,F,rho,mu); 
            break;
        case RNG_KE:   
            turb = new RNG_KE_Model(U,F,rho,mu); 
            break;
        case REALIZABLE_KE:   
            turb = new REALIZABLE_KE_Model(U,F,rho,mu); 
            break;
        case KW:   
            turb = new KW_Model(U,F,rho,mu); 
            break;
        case LES:  
            bneedWallDist = true;
            turb = new LES_Model(U,F,rho,mu); 
            break;
        default:
            turb = new Turbulence_Model(U,F,rho,mu); 
            break;
    }
    turb->enroll();
    return turb;
}
