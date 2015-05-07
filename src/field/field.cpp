#include "field.h"

using namespace std;

namespace Mesh {
	VectorVertexField vC(false);
	VectorFacetField  fC(false);
	VectorCellField   cC(false);
	VectorFacetField  fN(false);
	ScalarCellField   cV(false);
	ScalarFacetField  fI(false);
	ScalarCellField   yWall(false);
	IntVector         gFO;
	IntVector         gFN;
	IntVector  probeCells;
	Int   		gBCSfield;
	Int   		gBCSIfield;
	Int   		gBFSfield;
}
namespace Controls {
	Scheme convection_scheme = HYBRID;
	Int TVDbruner = 0;
	NonOrthoScheme nonortho_scheme = OVER_RELAXED;
	TimeScheme time_scheme = EULER;
	Scalar implicit_factor = 1;
	Int runge_kutta = 1;
	Scalar blend_factor = Scalar(0.2);
	Scalar tolerance = Scalar(1e-5f);
	Scalar dt = Scalar(.1);
	Scalar SOR_omega = Scalar(1.7);
	Solvers Solver = PCG; 
	Preconditioners Preconditioner = SSOR;
	State state = STEADY;
	Int max_iterations = 500;
	Int write_interval = 20;
	Int start_step = 0;
	Int end_step = 2;
	Int n_deferred = 0;
	Int save_average = 0;
	Int print_time = 0;
	CommMethod parallel_method = BLOCKED;
	Vector gravity = Vector(0,-9.81,0);
}
/*
 * Load mesh
 */
bool Mesh::LoadMesh(Int step,bool first, bool remove_empty) {
	if(MP::printOn)
		MP::printH("Reading mesh:\n");
	if(readMesh(step,first)) {
		MP::printH("\t%d vertices\n\t%d facets\n\t%d cells\n",
			gVertices.size(),gFacets.size(),gCells.size());
		/*initialize mesh*/
		addBoundaryCells();
		calcGeometry();
		DG::init_poly();
		/* remove empty faces*/
		if(remove_empty) {
			Boundaries::iterator it = gBoundaries.find("delete");
			if(it != gBoundaries.end()) {
				IntVector& fs = gBoundaries["delete"];
				removeBoundary(fs);
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
		/*geometric mesh fields*/
		remove_fields();
		initGeomMeshFields();
		if(MP::printOn) 
			cout << "--------------------------------------------\n";
		return true;
	}
	return false;
}
/*
 * Initialize geometric mesh fields
 */
void Mesh::initGeomMeshFields() {
	gFO = gFOC;
	gFN = gFNC;
	gBCSfield = gBCS * DG::NP;
	gBCSIfield = gBCSI * DG::NP;
	/* Allocate fields*/
	vC.deallocate(false);
	vC.allocate(gVertices);
	fC.deallocate(false);
	fC.allocate();
	cC.deallocate(false);
	cC.allocate();
	fN.deallocate(false);
	fN.allocate();
	cV.deallocate(false);
	cV.allocate();
	fI.deallocate(false);
	fI.allocate();
	/*expand*/
	forEach(_cC,i) {
		cC[i] = _cC[i];
		cV[i] = _cV[i];
	}
	forEach(_fC,i) {
		fC[i] = _fC[i];
		fN[i] = _fN[i];
	}
	if(DG::NPMAT) {
		DG::expand(cC);
		DG::expand(cV);
		DG::expand(fC);
		DG::expand(fN);
		DG::expand(fI);
		DG::init_basis();
	}
	/*Ghost face marker*/
	IntVector isGhostFace;
	isGhostFace.assign(gFacets.size(),0);
	forEach(gInterMesh,i) {
		interBoundary& b = gInterMesh[i];
		forEach(*b.f,j) {
			Int faceid = (*b.f)[j];
			isGhostFace[faceid] = 1;
		}
	}
	/* Facet interpolation factor to the owner of the face.
	 * Neighbor takes (1 - f) */
	forEach(gFacets,faceid) {
		Vector v0 = gVertices[gFacets[faceid][0]];
		for(Int n = 0; n < DG::NPF;n++) {
			Int k = faceid * DG::NPF + n;
			Int c1 = gFO[k];
			Int c2 = gFN[k];
			if(c2 >= gBCSfield && !isGhostFace[faceid]) {
				fI[k] = 0;
				cV[c2] = cV[c1];
				cC[c2] = fC[k];
			} else if(equal(cC[c1],fC[k]))  
				fI[k] = 0.5;
			else 
				fI[k] = 1 - dot(v0 - cC[c1],fN[k]) / 
					        dot(cC[c2] - cC[c1],fN[k]);
		}
	}
	/*exchange ghost cell volume and center*/
	exchange_ghost(&cV[0]);
	exchange_ghost(&cC[0]);
	/*Construct wall distance field*/
	{
		yWall.deallocate(false);
		yWall.construct("yWall");
		yWall = Scalar(0);
		//boundary
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
		applyExplicitBCs(yWall,true,true);
	}
}
/*find nearest cell*/
Int Mesh::findNearestCell(const Vector& v) {
	Scalar mindist,dist;
	Int bi = 0;
	mindist = mag(v - cC[0]);
	for(Int i = 0;i < gBCSfield;i++) {
		dist = mag(v - cC[i]);
		if(dist < mindist) {
			mindist = dist;
			bi = i;
		}
	}
	return bi;
}
Int Mesh::findNearestFace(const Vector& v) {
	Scalar mindist,dist;
	Int bi = 0;
	mindist = mag(v - fC[0]);
	forEach(fC,i) {
		dist = mag(v - fC[i]);
		if(dist < mindist) {
			mindist = dist;
			bi = i;
		}
	}
	return bi;
}
void Mesh::getProbeCells(IntVector& probes) {
	forEach(probePoints,j) {
		Vector v = probePoints[j];
		Int index = findNearestCell(v);
		probes.push_back(index);
	}
}
void Mesh::getProbeFaces(IntVector& probes) {
	forEach(probePoints,j) {
		Vector v = probePoints[j];
		Int index = findNearestFace(v);
		probes.push_back(index);
	}
}
void Mesh::calc_courant(const VectorCellField& U, Scalar dt) {
	ScalarCellField Courant;
	Courant = mag(U) * dt / pow(cV,1.0/3);
	Scalar minc = 1.0e200, maxc = 0.0;
	forEach(Courant,i) {
		if(Courant[i] < minc) minc = Courant[i];
		if(Courant[i] > maxc) maxc = Courant[i];
	}
	Scalar globalmax, globalmin;
	MP::allreduce(&maxc,&globalmax,1,MP::OP_MAX);
	MP::allreduce(&minc,&globalmin,1,MP::OP_MIN);
	if(MP::printOn) {
		MP::printH("Courant number: Max: %g Min: %g\n",
			globalmax,globalmin);
	}
}
/*
 * Read/Write
 */
void Mesh::write_fields(Int step) {
	forEachCellField(writeAll(step));
}
void Mesh::read_fields(Int step) {
	forEachCellField(readAll(step));
}
void Mesh::remove_fields() {
	forEachCellField(removeAll());
	forEachFacetField(removeAll());
	forEachVertexField(removeAll());
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
	params.enroll("implicit_factor",&implicit_factor);

	params.enroll("probe",&Mesh::probePoints);
	
	params.enroll("gravity", &gravity);
	
	Option* op;
	op = new Option(&convection_scheme,17,
		"CDS","UDS","HYBRID","BLENDED","LUD","CDSS","MUSCL","QUICK",
		"VANLEER","VANALBADA","MINMOD","SUPERBEE","SWEBY","QUICKL","UMIST",
		"DDS","FROMM");
	params.enroll("convection_scheme",op);
	op = new BoolOption(&TVDbruner);
	params.enroll("tvd_bruner",op);
	op = new Option(&nonortho_scheme,4,"NONE","MINIMUM","ORTHOGONAL","OVER_RELAXED");
	params.enroll("nonortho_scheme",op);
	op = new Option(&time_scheme,2,"EULER","SECOND_ORDER");
	params.enroll("time_scheme",op);
	params.enroll("runge_kutta",&runge_kutta);
	op = new Option(&Solver,3,"JACOBI","SOR","PCG");
	params.enroll("method",op);
	op = new Option(&Preconditioner,4,"NONE","DIAG","SSOR","DILU");
	params.enroll("preconditioner",op);
	op = new Option(&state,2,"STEADY","TRANSIENT");
	params.enroll("state",op);
	op = new Option(&parallel_method,2,"BLOCKED","ASYNCHRONOUS");
	params.enroll("parallel_method",op);
	op = new Util::BoolOption(&save_average);
	params.enroll("average",op);
	params.enroll("print_time",&print_time);
	params.enroll("npx",&DG::Nop[0]);
	params.enroll("npy",&DG::Nop[1]);
	params.enroll("npz",&DG::Nop[2]);
}

