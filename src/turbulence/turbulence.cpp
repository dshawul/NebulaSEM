#include "turbulence.h"
#include "ke.h"
#include "rngke.h"
#include "realizableke.h"
#include "kw.h"
#include "les.h"

Int Turbulence_Model::turb_model = 2;
bool Turbulence_Model::bneedWallDist = false;

void Turbulence_Model::RegisterTable(Util::ParamList& params) {
	Util::Option* op;
	op = new Util::Option(&turb_model,7,
		"NONE","MIXING_LENGTH","KE","RNG_KE","REALIZABLE_KE","KW","LES");
	params.enroll("turbulence_model",op);
}

Turbulence_Model* Turbulence_Model::Select(VectorCellField& U,ScalarFacetField& F,
					Scalar& rho,Scalar& viscosity,bool& Steady) {
	/*turbulence model*/
	enum TurbModel {
		NONE,MIXING_LENGTH,KE,RNG_KE,REALIZABLE_KE,KW,LES
	};
	bneedWallDist = false;
	
	/*Select turbulence model*/
	Turbulence_Model* turb;
	switch(turb_model) {
		case KE:   
			turb = new KE_Model(U,F,rho,viscosity,Steady); 
			break;
		case MIXING_LENGTH:   
			bneedWallDist = true;
			turb = new MixingLength_Model(U,F,rho,viscosity,Steady); 
			break;
		case RNG_KE:   
			turb = new RNG_KE_Model(U,F,rho,viscosity,Steady); 
			break;
		case REALIZABLE_KE:   
			turb = new REALIZABLE_KE_Model(U,F,rho,viscosity,Steady); 
			break;
		case KW:   
			turb = new KW_Model(U,F,rho,viscosity,Steady); 
			break;
		case LES:  
			bneedWallDist = true;
			turb = new LES_Model(U,F,rho,viscosity,Steady); 
			break;
		default:
			turb = new Turbulence_Model(U,F,rho,viscosity,Steady); 
			break;
	}
	turb->enroll();
	return turb;
}