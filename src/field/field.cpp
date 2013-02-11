#include "field.h"

using namespace std;

namespace Mesh {
	VectorVertexField vC;
	VectorFacetField  fC;
	VectorCellField   cC;
	VectorFacetField  fN;
	ScalarCellField   cV;
	ScalarFacetField  fI;
	ScalarCellField   yWall(false);
}
namespace Controls {
	Scheme convection_scheme = HYBRID;
	Int TVDbruner = 0;
	Scheme interpolation_scheme = CDS;
	NonOrthoScheme nonortho_scheme = OVER_RELAXED;
	Scalar time_scheme_factor = 1;
	Scalar blend_factor = Scalar(0.2);
	Scalar tolerance = Scalar(1e-5f);
	Scalar dt = Scalar(.1);
	Scalar SOR_omega = Scalar(1.7);
	Solvers Solver = PCG; 
	Preconditioners Preconditioner = SORP;
	State state = STEADY;
	Int max_iterations = 500;
	Int write_interval = 20;
	Int start_step = 0;
	Int end_step = 2;
	Int n_deferred = 0;
	Int save_average = 0;
	CommMethod ghost_exchange = BLOCKED;
	CommMethod parallel_method = BLOCKED;
}

/*
 * Initialize geometric mesh fields
 */
void Mesh::initGeomMeshFields(bool remove_empty) {
	/*initialize mesh*/
	addBoundaryCells();
	calcGeometry();
	/* remove empty faces*/
	if(remove_empty) {
		Boundaries::iterator it = gBoundaries.find("delete");
		if(it != gBoundaries.end()) {
			removeBoundary(gBoundaries["delete"]);
			gBoundaries.erase(it);
		}
	}
	/*erase interior and empty boundaries*/
	for(Boundaries::iterator it = gBoundaries.begin();
                it != gBoundaries.end();) {
		if(it->second.size() <= 0 || 
			it->first.find("interior") != std::string::npos
			) {
				gBoundaries.erase(it++);
		} else ++it;
	}
	/* Allocate fields*/
	vC.allocate(gVertices);
	fC.allocate(_fC);
	cC.allocate(_cC);
	fN.allocate(_fN);
	cV.allocate(_cV);
	fI.allocate();
	/* Facet interpolation factor to the owner of the face.
	 * Neighbor takes (1 - f) */
	exchange_ghost(&cV[0]);
	exchange_ghost(&cC[0]);
	forEach(gFacets,i) {
		Int c1 = gFO[i];
		Int c2 = gFN[i];
		Scalar s1 = mag(cC[c1] - fC[i]);
		Scalar s2 = mag(cC[c2] - fC[i]);
		fI[i] = 1.f - s1 / (s1 + s2);
	}
	/*Construct wall distance field*/
	{
		yWall.construct("yWall");
		yWall = Scalar(0);
		/*boundary*/
		BCondition<Scalar>* bc;
		forEachIt(Boundaries,gBoundaries,it) {
			string bname = it->first;
			bc = new BCondition<Scalar>(yWall.fName);
			bc->bname = bname;
			if(bname.find("WALL") != std::string::npos) {
				bc->cname = "DIRICHLET";
				bc->value = Scalar(0);
			} else if(bname.find("interMesh") != std::string::npos) {
			} else {
				bc->cname = "NEUMANN";
				bc->value = Scalar(0);
			}
			bc->init_indices();
			AllBConditions.push_back(bc);
		}
		updateExplicitBCs(yWall,true,true);
	}
}
/*
 * Read/Write
 */
void Mesh::write_fields(Int step) {
	forEachField(writeAll(step));
}
void Mesh::read_fields(Int step) {
	forEachField(readAll(step));
}
void Mesh::enroll(Util::ParamList& params) {
	using namespace Controls;
	using namespace Util;

	params.enroll("max_iterations",&max_iterations);
	params.enroll("write_interval",&write_interval);
	params.enroll("start_step",&start_step);
	params.enroll("end_step",&end_step);
	params.enroll("n_deferred",&n_deferred);

    params.enroll("blend_factor",&blend_factor);
	params.enroll("tolerance",&tolerance);
	params.enroll("dt",&dt);
	params.enroll("SOR_omega",&SOR_omega);
	params.enroll("time_scheme_factor",&time_scheme_factor);

	params.enroll("probe",&Mesh::probePoints);

	Option* op;
	op = new Option(&convection_scheme,17,
		"CDS","UDS","HYBRID","BLENDED","LUD","CDSS","MUSCL","QUICK",
		"VANLEER","VANALBADA","MINMOD","SUPERBEE","SWEBY","QUICKL","UMIST",
		"DDS","FROMM");
	params.enroll("convection_scheme",op);
	op = new BoolOption(&TVDbruner);
	params.enroll("tvd_bruner",op);
	op = new Option(&interpolation_scheme,2,"CDS","UDS");
	params.enroll("interpolation_scheme",op);
	op = new Option(&nonortho_scheme,4,"NONE","MINIMUM","ORTHOGONAL","OVER_RELAXED");
	params.enroll("nonortho_scheme",op);
	op = new Option(&Solver,3,"JACOBI","SOR","PCG");
	params.enroll("method",op);
	op = new Option(&Preconditioner,4,"NONE","DIAG","SOR","DILU");
	params.enroll("preconditioner",op);
	op = new Option(&state,2,"STEADY","TRANSIENT");
	params.enroll("state",op);
	op = new Option(&ghost_exchange,2,"BLOCKED","ASYNCHRONOUS");
	params.enroll("ghost_exchange",op);
	op = new Option(&parallel_method,2,"BLOCKED","ASYNCHRONOUS");
	params.enroll("parallel_method",op);
	op = new Util::BoolOption(&save_average);
	params.enroll("average",op);
}
