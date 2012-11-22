#include "field.h"

using namespace std;

namespace Mesh {
	VectorVertexField vC;
	VectorFacetField  fC;
	VectorCellField   cC;
	VectorFacetField  fN;
	ScalarCellField   cV;
	ScalarFacetField  fI;
}
namespace Controls {
	Scheme convection_scheme = HYBRID;
	HigherOrderScheme higher_scheme = IMPLICIT;
	Scheme interpolation_scheme = CDS;
	NonOrthoScheme nonortho_scheme = OVER_RELAXED;
	Scalar time_scheme_factor = 1;
	Scalar blend_factor = Scalar(0.2);
	Scalar tolerance = Scalar(1e-5f);
	Scalar dt = Scalar(.1);
	Scalar SOR_omega = Scalar(1.7);
	Solvers Solver = PCG; 
	State state = STEADY;
	Int max_iterations = 500;
	Int write_interval = 20;
	Int start_step = 0;
	Int end_step = 2;
	GhostExchange ghost_exchange = BLOCKED;
}
namespace {
	vector<Vector> _fC;
	vector<Vector> _cC;
	vector<Vector> _fN;
	vector<Scalar> _cV;
	vector<Scalar> _fI;
}
/* 
* Remove empty boundary
*/
void Mesh::removeBoundary(IntVector& fs) {
	cout << "Removing faces: " << fs.size() << endl;

	Int i,j,count;
	IntVector Idf(gFacets.size(),0);
	IntVector Idc(gCells.size(),0);

	/*erase facet reference*/
	for(i = 0;i < fs.size();i++) {
		Int f = fs[i];
		Cell& co = gCells[gFO[f]];
		for(j = 0;j < co.size();j++) {
			if(co[j] == f) {
				co.erase(co.begin() + j); 
				break; 
			}
		}
		Cell& cn = gCells[gFN[f]];
		for(j = 0;j < cn.size();j++) {
			if(cn[j] == f) { 
				cn.erase(cn.begin() + j); 
				break; 
			}
		}
	}
	/*updated facet id*/
	for(i = 0;i < fs.size();i++)
		Idf[fs[i]] = Constants::MAX_INT;
	count = 0;
	for(i = 0;i < gFacets.size();i++) {
		if(Idf[i] != Constants::MAX_INT) 
			Idf[i] = count++;
		else
			gFacets[i].clear();
	}
	/*erase facets*/
	for(i = 0;i < gFacets.size();i++) {
		if(gFacets[i].size() == 0) {
			gFacets.erase(gFacets.begin() + i);
			gFO.erase(gFO.begin() + i);
			gFN.erase(gFN.begin() + i);
			_fC.erase(_fC.begin() + i);
			_fN.erase(_fN.begin() + i);
			_fI.erase(_fI.begin() + i);
			--i;
		}
	}
	/*updated facet id*/
	count = 0;
	for(i = 0;i < gCells.size();i++) {
		if(gCells[i].size() != 0) 
			Idc[i] = count++;
		else
			Idc[i] = Constants::MAX_INT;
	}
	/*erase cells*/
	for(i = 0;i < gCells.size();i++) {
		if(gCells[i].size() == 0) { 
			gCells.erase(gCells.begin() + i); 
			_cC.erase(_cC.begin() + i);
			_cV.erase(_cV.begin() + i);
			--i;
		} else {
			for(j = 0;j < gCells[i].size();j++) {
				gCells[i][j] = Idf[gCells[i][j]];
			}
		}
	}
	/*facet owner and neighbor*/
	for(i = 0;i < gFacets.size();i++) {
		gFO[i] = Idc[gFO[i]];
		gFN[i] = Idc[gFN[i]];
	}
	/*patches*/
	for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it) {
		IntVector& gB = it->second;
		for(i = 0;i < gB.size();i++) 
			gB[i] = Idf[gB[i]];
	}

	cout << "Total faces: " << gFacets.size() << endl;
}
/*
 * Initialize geometric mesh fields
 */
void Mesh::initGeomMeshFields() {
	Int i,j;

	/*first add boundary cells*/
	Mesh::addBoundaryCells();

	/*allocate*/
	_fC.assign(gFacets.size(),Vector(0));
	_cC.assign(gCells.size(),Vector(0));
	_fN.assign(gFacets.size(),Vector(0));
	_cV.assign(gCells.size(),Scalar(0));
	_fI.assign(gFacets.size(),Scalar(0));

	/* face centre*/
	for(i = 0;i < gFacets.size();i++) {
		Facet& f = gFacets[i];
		Vector C(0);
		for(j = 0;j < f.size();j++)
			C += gVertices[f[j]];
		_fC[i] = C / Scalar(f.size());
	}

	/* cell centre */
	for(i = 0;i < gCells.size();i++) {
		Cell& c = gCells[i];
		Vector C(0);
		for(j = 0;j < c.size();j++)
			C += _fC[c[j]];
		_cC[i] = C / Scalar(c.size());
	}
	/* face normal */
	Vector v1,v2,v3,v;
	for(i = 0;i < gFacets.size();i++) {
		Facet& f = gFacets[i];
		Vector N(0),C(0),Ni;
		v1 = _fC[i];
		for(j = 0;j < f.size();j++) {
			v2 = gVertices[f[j]];
			if(j + 1 == f.size())
				v3 = gVertices[f[0]];
			else
				v3 = gVertices[f[j + 1]];
			Ni = ((v2 - v1) ^ (v3 - v1));
			C += mag(Ni) * (v1 + v2 + v3) / 3;
			N += Ni;
		}
		_fC[i] = C / mag(N);    /*corrected face centre*/
		v = _fC[i] - _cC[gFO[i]];
		if((v & N) < 0) N = -N;
		_fN[i] = N / Scalar(2);
	}
	/* cell volumes */
	for(i = 0;i < gBCellsStart;i++) {
		Cell& c = gCells[i];
		Scalar V(0),Vi;
		Vector v = _cC[i],C(0);
		for(j = 0;j < c.size();j++) {
			v = _cC[i] - _fC[c[j]];
			Vi = mag(v & _fN[c[j]]);
			C += Vi * (2 * _fC[c[j]] + _cC[i]) / 3;
			V += Vi;
		}
		_cC[i] = C / V;         /*corrected cell centre */
		_cV[i] = V / Scalar(3);
	}
	for(i = gBCellsStart;i < gCells.size();i++) {
		_cV[i] = _cV[gFO[gCells[i][0]]];
	}
	/* Exchange ghost cells cC and cV */
	exchange_ghost(&_cV[0]);
	exchange_ghost(&_cC[0]);
	/* Facet interpolation factor to the owner of the face.
	 * Neighbor takes (1 - f) */
	for(i = 0;i < gFacets.size();i++) {
		Int c1 = gFO[i];
		Int c2 = gFN[i];
		Scalar s1 = mag(_cC[c1] - _fC[i]);
		Scalar s2 = mag(_cC[c2] - _fC[i]);
		_fI[i] = 1.f - s1 / (s1 + s2);
	}
	/* remove empty faces*/
	if(gBoundaries.find("delete") != gBoundaries.end())
		removeBoundary(gBoundaries["delete"]);
	/*
	 * Allocate fields
	 */
	vC.allocate(gVertices);
	fC.allocate(_fC);
	cC.allocate(_cC);
	fN.allocate(_fN);
	cV.allocate(_cV);
	fI.allocate(_fI);
}
/*
 * Read/Write
 */
void Mesh::write_fields(Int step) {
	ScalarCellField::write(step);
	VectorCellField::write(step);
	STensorCellField::write(step);
	TensorCellField::write(step);
}
void Mesh::read_fields(Int step) {
	ScalarCellField::read(step);
	VectorCellField::read(step);
	STensorCellField::read(step);
	TensorCellField::read(step);
}
void Mesh::enroll() {
	using namespace Controls;
	using namespace Util;

	StringParams::enroll("mesh",&gMeshName);

	IntParams::enroll("max_iterations",&max_iterations);
	IntParams::enroll("write_interval",&write_interval);
	IntParams::enroll("start_step",&start_step);
	IntParams::enroll("end_step",&end_step);

    ScalarParams::enroll("blend_factor",&blend_factor);
	ScalarParams::enroll("tolerance",&tolerance);
	ScalarParams::enroll("dt",&dt);
	ScalarParams::enroll("SOR_omega",&SOR_omega);
	ScalarParams::enroll("time_scheme_factor",&time_scheme_factor);

	VerticesParams::enroll("probe",&Mesh::probePoints);

	Option* op;
	op = new Option(&convection_scheme,14,
		"CDS","UDS","HYBRID","BLENDED","LUD","MUSCL","QUICK",
		"VANLEER","VANALBADA","MINMOD","SUPERBEE","SWEBY","QUICKL","UMIST");
	OptionParams::enroll("convection_scheme",op);
	op = new Option(&higher_scheme,2,"IMPLICIT","DEFERRED");
	OptionParams::enroll("higher_scheme",op);
	op = new Option(&interpolation_scheme,2,"CDS","UDS");
	OptionParams::enroll("interpolation_scheme",op);
	op = new Option(&nonortho_scheme,4,"NONE","MINIMUM","ORTHOGONAL","OVER_RELAXED");
	OptionParams::enroll("nonortho_scheme",op);
	op = new Option(&Solver,2,"SOR","PCG");
	OptionParams::enroll("method",op);
	op = new Option(&state,2,"STEADY","TRANSIENT");
	OptionParams::enroll("state",op);
	op = new Option(&ghost_exchange,2,"BLOCKED","ASYNCHRONOUS");
	OptionParams::enroll("ghost_exchange",op);
}
