#ifndef __TURBULENCE_H
#define __TURBULENCE_H

#include "field.h"
#include "solve.h"
/**
     \verbatim
     Description of RANS turbulence models
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Navier Stokes without source terms:
         d(rho*u)/dt + div(rho*uu) = -grad(p) + div(V)
           where u,rho are instantaneous quantities, and V is the 
           deviatoric stress tensor (a.k,a Viscous stress tensor ).
     RANS (Reynolds averaged navier stokes): with mean quantites U, rho
         d(rho*U)/dt + div(rho*UU) + div(rho*u'u') = -grad(P) + div(V)
         d(rho*U)/dt + div(rho*UU) = -grad(P) + div(V + R)
     where the Viscous (V) and Reynolds (R) stress tensors are
         V =  2*mu*dev(S)
         R = -rho*u'u'
       where S = (gU + gUt) / 2 is the mean strain tensor.
         2 * dev(S) = (gU + gUt) - 2 * gUt_ii / 3
                    = gU + (gUt - 2 * gUt_ii/3)
                    = gU + dev(gUt,2)
     The viscouse and reynolds stresses contain normal components that could
     be absored to the pressure term for computational efficiency reasons.
              V = 2 * mu * dev(S)
     Boussinesq model for the deviatoric Reynolds stress uses eddy viscosity (emu):
         dev(R) = 2 * emu * dev(S) 
              R = 2 * emu * dev(S) + R_ii/3
     Viscous and Reynolds stress together:
         V + R = 2 * (mu + emu) * dev(S) + R_ii/3
               = 2 * eff_mu * dev(S) + R_ii/3
     Volume integrated V+R (i.e. force), results in a laplacian for the implict term
         div(V + R) = div(2 * eff_mu * dev(S)) + div(R_ii/3)
                    = div(eff_mu*gU) + div(eff_mu * dev(gUt,2))  + div(R_ii/3)
                         Implicit           Explicit               Absorbed in pressure
                                                                    P_m = P - R_ii/3
     For constant effective viscosity (i.e. no turbulence and constant dynamic viscosity)
              (Implicit) div(eff_mu * gU)  = eff_mu * div(grad(U)) = eff_mu * lapalacian(U)
              (Explicit) div(eff_mu * gUt) = eff_mu * grad(div(U)).
     In this case, the explicit term can be dropped for incompressible flows where div(U) = 0.

     Final RANS equation after substituting div(V+R):
         d(rho*U)/dt + div(rho*UU) = -grad(P) + div(V + R)
                                   = -grad(P_m) + div(eff_mu*gU) + div(eff_mu * dev(gUt,2))
     where the modified pressure P_m, a.k.a mechanical pressure, can be calculated from 
     thermodynamic pressure P as
         P_m = P - R_ii/3
             = P + 2/3*rho*k
     since R_ii/3 = trace(-rho*u'u') = -2/3*rho*k*I

     Base turbulence model:
         This default class has no turbulence model so it is a laminar solver. 
         Only the viscous stress V is added to the NS equations. Turbulence models 
         derived from this class add a model for Reynold's stress R usually by solving 
         some turbulence transport equations.
     \endverbatim
*/
struct Turbulence_Model {

    VectorCellField& U;
    VectorCellField& Fc;
    ScalarFacetField& F;
    ScalarCellField& rho;
    ScalarCellField& mu;

    Util::ParamList params;
    bool writeStress;
    /*constructor*/
    Turbulence_Model(VectorCellField& tU,VectorCellField& tFc,ScalarFacetField& tF,ScalarCellField& trho,ScalarCellField& tmu) :
        U(tU),
        Fc(tFc),
        F(tF),
        rho(trho),
        mu(tmu),
        params("turbulence"),
        writeStress(false)
    {
    }
    virtual ~Turbulence_Model() {};
    /*overridable functions*/
    virtual void enroll() {
        using namespace Util;
        Option* op = new BoolOption(&writeStress);
        params.enroll("writeStress",op);
    };
    virtual void solve() {};
    /* Turbulence viscosity */
    virtual ScalarCellField getTurbVisc() {
        return Scalar(0);
    }
    /* Explicit V */
    virtual VectorCellField getExplicitStresses() {
        return divi(mu * dev(trn(gradi(U)),2.0));
    };
    /* V */
    STensorCellField getViscousStress() {
        return 2.0 * mu * dev(sym(gradi(U)),1.0);
    }
    /* R */
    virtual STensorCellField getReynoldsStress() {
        return STensor(0);
    }
    /* TKE */
    virtual ScalarCellField getK() {
        return Scalar(0);
    }
    /* Turbulence model selection */
    static Int turb_model;
    static bool bneedWallDist;
    static bool needWallDist() { return bneedWallDist;}
    static void RegisterTable(Util::ParamList& params);
    static Turbulence_Model* Select(VectorCellField& U,VectorCellField& Fc,ScalarFacetField& F,
        ScalarCellField& rho,ScalarCellField& mu);
};
/**
 * Eddy viscosity models based on Boussinesq's assumption
 * that the action of Reynolds and Viscous stress are similar.
 */
struct EddyViscosity_Model : public Turbulence_Model {
    ScalarCellField eddy_mu; 
    enum Model {
        SMAGORNSKY,BALDWIN,KATO
    };
    enum WallModel {
        NONE,STANDARD,LAUNDER
    };
    Model modelType;
    WallModel wallModel;
    
    /*constructor*/
    EddyViscosity_Model(VectorCellField& tU,VectorCellField& tFc,ScalarFacetField& tF,ScalarCellField& trho,ScalarCellField& tmu) :
        Turbulence_Model(tU,tFc,tF,trho,tmu),
        eddy_mu("emu",READWRITE),
        modelType(SMAGORNSKY),
        wallModel(STANDARD)
    {
    }
    /*Register options*/
    virtual void enroll() {
        using namespace Util;
        Option* op = new Option(&modelType,
            {"SMAGORNSKY","BALDWIN","KATO"});
        params.enroll("modelType",op);
        Turbulence_Model::enroll();
    }
    /*eddy_mu*/
    virtual void calcEddyViscosity(const TensorCellField& gradU) = 0;
    virtual ScalarCellField getTurbVisc() {
        return eddy_mu;
    }
    /* Explicit V + R */
    virtual VectorCellField getExplicitStresses() {
        TensorCellField gradU = gradi(U);
        calcEddyViscosity(gradU);
        setWallEddyMu();
        fillBCs(eddy_mu);
        return divi((eddy_mu + mu) * dev(trn(gradU),2.0));
    };

    /* R */
    virtual STensorCellField getReynoldsStress() {
        STensorCellField R = 2 * eddy_mu * dev(sym(gradi(U)),1.0) - 
                 STensorCellField(Constants::I_ST) * (2 * rho * getK() / 3);
        return R;
    }
    /* S2 */
    ScalarCellField getS2(const TensorCellField& gradU) {
        ScalarCellField magS;
        if(modelType == SMAGORNSKY) {
            STensorCellField S = sym(gradU);
            magS = S & S;
        } else if(modelType == BALDWIN) {
            TensorCellField O = skw(gradU);
            magS = O & O;
        } else {
            STensorCellField S = sym(gradU);
            TensorCellField O = skw(gradU);
            magS = sqrt((S & S) * (O & O));
        }
        return (2 * magS);
    }
    /*Fix near wall cell values*/
    void FixNearWallValues(ScalarCellMatrix& M) {
        using namespace Mesh;
        BasicBCondition* bbc;
        forEach(AllBConditions,d) {
            bbc = AllBConditions[d];
            if(bbc->fIndex == eddy_mu.fIndex && bbc->cIndex == Mesh::ROUGHWALL) {
                IntVector& wall_faces = *bbc->bdry;
                if(wall_faces.size()) {
                    Int f,c1;
                    forEach(wall_faces,i) {
                        f = wall_faces[i];
                        c1 = FO[f];
                        M.Fix(c1,(*M.cF)[c1]);
                    }
                }
            }
        }
    }
    /* Wall functions */
    void setWallEddyMu() {
        using namespace Mesh;
        BasicBCondition* bbc;
        forEach(AllBConditions,d) {
            bbc = AllBConditions[d];
            if(bbc->fIndex == eddy_mu.fIndex && bbc->cIndex == Mesh::ROUGHWALL) {
                IntVector& wall_faces = *bbc->bdry;
                LawOfWall& low = bbc->low;
                if(wall_faces.size()) {
                    forEach(wall_faces,i) {
                        applyWallFunction(wall_faces[i],low);
                    }
                }
            }
        }
    }
    /*overridable*/
    virtual void applyWallFunction(Int f,LawOfWall& low) = 0;
};
/**
 * Base two equation K-X turbulence model
 */ 
struct KX_Model : public EddyViscosity_Model {
    /*model coefficients*/
    Scalar Cmu;
    Scalar SigmaK;
    Scalar SigmaX;
    Scalar C1x;
    Scalar C2x;

    Scalar k_UR;
    Scalar x_UR;

    /*turbulence fields*/
    ScalarCellField k;         
    ScalarCellField x;        
    ScalarCellField Pk;       

    /*constructor*/
    KX_Model(VectorCellField& tU,VectorCellField& tFc,ScalarFacetField& tF,ScalarCellField& trho,ScalarCellField& tmu,const char* xname) :
        EddyViscosity_Model(tU,tFc,tF,trho,tmu),
        k_UR(0.7),
        x_UR(0.7),
        k("k",READWRITE),
        x(xname,READWRITE)
    {
        wallModel = LAUNDER;
    }
    /*TKE*/
    virtual ScalarCellField getK() { return k; }
    /*Register options*/
    virtual void enroll() {
        using namespace Util;
        params.enroll("k_UR",&k_UR);
        params.enroll("x_UR",&x_UR);
        EddyViscosity_Model::enroll();
    }
    /* k-x model specific over-ridables*/
    virtual void calcEddyMu() = 0;
    virtual Scalar calcX(Scalar ustar,Scalar kappa,Scalar y) = 0;
    virtual Scalar getCmu(Int i) { 
        return Cmu; 
    }
    /* eddy viscosity*/
    virtual void calcEddyViscosity(const TensorCellField& gradU) {
        calcEddyMu();
        Pk = getS2(gradU) * eddy_mu;
    }
    /* wall function */
    virtual void applyWallFunction(Int f,LawOfWall& low) { 
        using namespace Mesh;
        Int c1 = FO[f];
        Int c2 = FN[f];

        /*calc ustar*/
        Scalar nu = mu[c1] / rho[c1];
        Scalar ustar = Scalar(0.0);
        Scalar y = mag(unit(fN[f]) & (cC[c1] - cC[c2]));
        if(wallModel == STANDARD) {
            ustar = low.getUstar(nu,mag(U[c1]),y);
            k[c1] = pow(ustar,2) / sqrt(getCmu(c1));
        } else if(wallModel == LAUNDER) {
            ustar = pow(getCmu(c1),Scalar(0.25)) * sqrt(k[c1]);
        }
        x[c1] = calcX(ustar,low.kappa,y);

        /* calculate eddy viscosity*/
        Scalar yp = (ustar * y) / nu;
        Scalar up = low.getUp(ustar,nu,yp);                                         
        eddy_mu[c1] = mu[c1] * (yp / up - 1);

        /* turbulence generation and dissipation */
        if(wallModel == LAUNDER) {
            Scalar mag_dudy = mag((U[c2] - U[c1]) / y);
            Scalar mag_dudy_log = ustar / (low.kappa * y);
            Pk[c1] = (mag_dudy * mag_dudy_log) * eddy_mu[c1];
        }
    };
};

#endif
