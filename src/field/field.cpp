#include "field.h"

using namespace std;

namespace Mesh {
	VectorVertexField vC(false);
	VectorFacetField  fC(false);
	VectorCellField   cC(false);
	VectorFacetField  fN(false);
	ScalarCellField   cV(false);
	ScalarFacetField  fI(false);
	ScalarFacetField  fD(false);
	ScalarCellField   yWall(false);
	IntVector         gFO;
	IntVector         gFN;
	IntVector  probeCells;
	Int   		gBCSfield;
	Int   		gBCSIfield;
	Int   		gBFSfield;
}

std::list<BaseField*> BaseField::allFields;
std::vector<std::string> BaseField::fieldNames;

namespace Controls {
	Scheme convection_scheme = HYBRID;
	Int TVDbruner = 0;
	NonOrthoScheme nonortho_scheme = OVER_RELAXED;
	TimeScheme time_scheme = BDF1;
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
	Int amr_step = 0;
	Int n_deferred = 0;
	Int save_average = 0;
	Int print_time = 0;
	CommMethod parallel_method = BLOCKED;
	Vector gravity = Vector(0,0,-9.860616);
	RefineParams refine_params;
}
/*find the last grid loaded*/
static int findLastRefinedGrid(Int step_) {
	int step = step_;
	for(;step >= 0;--step) {
		stringstream path;
		path << Mesh::gMeshName << "_" << step;
		string str = path.str();
		ifstream is(str.c_str());
		if(!is.fail())
			break;
	}
	if(step < 0) step = step_;
	return step;
}
/*
 * Load mesh
 */
bool Mesh::LoadMesh(Int step_,bool first, bool remove_empty, bool coarse) {
	/*load refined mesh*/
	int step = step_;
	if(first && !coarse)
		step = findLastRefinedGrid(step_);
	
	/*load mesh*/
	if(readMesh(step,first,coarse)) {
		if(MP::host_id == 0)
			cout << "--------------------------------------------\n";
		MP::printH("\t%d vertices\t%d facets\t%d cells\n",
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
		if(MP::host_id == 0) 
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
	fD.deallocate(false);
	fD.allocate();
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
	/*Start communicating cV and cC*/
	ASYNC_COMM<Scalar> commv(&cV[0]);
	ASYNC_COMM<Vector> commc(&cC[0]);
	commv.send();
	commc.send();
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
	/*finish comm*/
	commv.recv();
	commc.recv();
	/*construct diffusivity factor*/
	if(DG::NPMAT) {
		using namespace DG;
		
		//penalty
		Scalar k = (NPX > NPY) ? NPX : 
				  ((NPY > NPZ) ? NPY : NPZ);
		Scalar num = (k + 1) * (k + 3) / 3;
		
		fD = cds(cV);
		forEach(fN,i)
			fD[i] = num * mag(fN[i]) / (fD[i]);
		
		//diffusivity
		VectorCellField grad_psi = Vector(0);

		for(Int ci = 0; ci < gBCS; ci++) {
			forEachLglBound(ii,jj,kk) {
				Int index = INDEX4(ci,ii,jj,kk);
				
#define PSID(im,jm,km) {					\
	Int index1 = INDEX4(ci,im,jm,km);		\
	Vector dpsi_ij;							\
	DPSIR(dpsi_ij,im,jm,km);				\
	dpsi_ij = dot(dpsi_ij,Jinv[index1]);	\
	grad_psi[index] += dpsi_ij;				\
}
				forEachLglX(i) PSID(i,jj,kk);
				forEachLglY(j) if(j != jj) PSID(ii,j,kk);
				forEachLglZ(k) if(k != kk) PSID(ii,jj,k);
			}
		}

		fD += dot(cds(grad_psi),fN);
			
	} else {
		using namespace Controls;
		
		forEach(fD,i) {
			Int c1 = gFO[i];
			Int c2 = gFN[i];
			Vector dv = cC[c2] - cC[c1];
			if(nonortho_scheme == OVER_RELAXED) {
				fD[i] = ((fN[i] & fN[i]) / (fN[i] & dv));
			} else if(nonortho_scheme == MINIMUM) {
				fD[i] = ((fN[i] & dv) / (dv & dv));
			} else {
				fD[i] = sqrt((fN[i] & fN[i]) / (dv & dv));
			}
		}
	}
	/*Construct wall distance field*/
	{
		yWall.deallocate(false);
		yWall.construct();
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
	BaseField::allFields.clear();
}
void Mesh::enroll(Util::ParamList& params) {
	using namespace Controls;
	using namespace Util;

	params.enroll("max_iterations",&max_iterations);
	params.enroll("write_interval",&write_interval);
	params.enroll("start_step",&start_step);
	params.enroll("end_step",&end_step);
	params.enroll("amr_step",&amr_step);
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
	op = new Option(&time_scheme,6,"BDF1","BDF2","BDF3","BDF4","BDF5","BDF6");
	params.enroll("time_scheme",op);
	params.enroll("runge_kutta",&runge_kutta);
	op = new Option(&Solver,3,"JAC","SOR","PCG");
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
void Controls::enrollRefine(Util::ParamList& params) {
	Util::Option* op;
	op = new Util::Option(&refine_params.shape, 2,"QUAD","TRI");
	params.enroll("shape",op);
	params.enroll("direction",&refine_params.dir);
	params.enroll("field",&refine_params.field);
	params.enroll("field_max",&refine_params.field_max);
	params.enroll("field_min",&refine_params.field_min);
	params.enroll("limit",&refine_params.limit);
}
/*
* Adaptive mesh refinement
*/
void Prepare::refineMesh(const RefineParams& rparams, Int step) {
	using namespace Mesh;
	
	/*create fields*/
	Prepare::createFields(BaseField::fieldNames,step);
	Prepare::readFields(BaseField::fieldNames,step);
	
	/*Unrefine fields and mesh*/
	{
		IntVector cellMap;
		stringstream path;
		int stepn = findLastRefinedGrid(step);
		path << "cellMap" << "_" << stepn;
		ifstream is(path.str().c_str());
		if(!is.fail()) {
			is >> hex;
			is >> cellMap;
			is >> dec;

			forEachIt(std::list<BaseField*>, BaseField::allFields, it)
				(*it)->unrefineField(step,cellMap);
		}
		
		/*load coarse mesh*/
		LoadMesh(step,true,false,true);
	}
	
	/*create fields*/
	Prepare::createFields(BaseField::fieldNames,step);
	Prepare::readFields(BaseField::fieldNames,step);
	
	/*find cells to refine*/
	IntVector rCells;
	{
		/*find quantity of interest*/
		ScalarCellField qoi;
		BaseField* bf = BaseField::findField(rparams.field);
		if(bf) bf->norm(&qoi);
		// qoi = mag(gradi(qoi));
	
		/*max and min*/
		Scalar maxq = 0,minq = 10e30;
		for(Int i = 0;i < gBCS;i++) {
			if(qoi[i] > maxq) maxq = qoi[i];
			if(qoi[i] < minq) minq = qoi[i];
		}
		maxq = max(Constants::MachineEpsilon,maxq);
		qoi /= maxq;
		maxq = 1;
		minq = minq / maxq;
	
		/*get cells to refine*/
		gCells.erase(gCells.begin() + gBCS,gCells.end());
		for(Int i = 0;i < gBCS;i++) {
			if(qoi[i] >= rparams.field_max)
				rCells.push_back(i);
		}
	}
	
	/*refine mesh and fields*/
	{
		IntVector cellMap;
		refineMesh(rCells,rparams,cellMap);
		forEachIt(std::list<BaseField*>, BaseField::allFields, it)
			(*it)->refineField(step,cellMap); 

		/*Write cellMap*/
		{
			stringstream path;
			path << "cellMap" << "_" << step;
			ofstream os(path.str().c_str());
			os << hex;
			os << cellMap;
			os << dec;
		}

		/*Write mesh*/
		{
			stringstream path;
			path << gMeshName << "_" << step;
			ofstream os(path.str().c_str());
			gMesh.write(os);
		}
	}
	
	/*destroy*/ 
	BaseField::destroyFields();
}

void Prepare::refineMesh(IntVector& rCells,const RefineParams& rparams,IntVector& cellMap) {
	using namespace Mesh;
			
	/*build refine flags array*/
	IntVector refineC,refineF,
			  rFacets,rsFacets;
	cellMap.assign(gBCS,0);
	forEach(cellMap,i)
		cellMap[i] = i;
	refineC.assign(gBCS,0);
	refineF.assign(gFacets.size(),0);
	forEach(rCells,i) {
		Int ci = rCells[i];
		refineC[ci] = 1;
		Cell& c = gCells[ci];
		forEach(c,j)
			refineF[c[j]] = 1;
	}
	
	/*choose faces to refine*/
	bool refine3D = equal(Scalar(0.0),mag(rparams.dir));
	forEach(refineF,i) {
		if(refineF[i]) {
			Vector eu = unit(fN[i]);
			if(refine3D || 
				equal(rparams.dir,eu) || 
				equal(rparams.dir,-eu))
				rFacets.push_back(i);
			else
				refineF[i] = 0;
		}
	}
	
	/***************************
	 * Refine facets
	 ***************************/
	Int nj;
	IntVector startF;
	startF.assign(gFacets.size(),0);
	Int ivBegin = gVertices.size();
	forEach(rFacets,i) {
		Int fi = rFacets[i];
		Facet f = gFacets[fi];
		startF[fi] = gFacets.size();
		Vector C = fC[fi];
		gVertices.push_back(C);
		Int fci = gVertices.size() - 1;
		
		/*quadrilaterals*/
		if(rparams.shape == 0) {
			/*add vertices*/
			IntVector midpts;
			Vector v1,v2;
			forEach(f,j) {
				v1 = gVertices[f[j]];
				nj = j + 1;
				if(nj == f.size())
					nj = 0;
				v2 = gVertices[f[nj]];
				Vector Ce = (v1 + v2) / 2.0;
			
				//duplicate
				Int k = ivBegin;
				for(; k < gVertices.size();k++) {
					if(equal(Ce,gVertices[k])) {
						midpts.push_back(k);
						break;
					}
				}
				if(k == gVertices.size()) {
					gVertices.push_back(Ce);
					midpts.push_back(k);
				}
			}
		
			/*add facets*/
			Int k1,k2;
			forEach(f,j) {
				k1 = midpts[j];
				if(j == 0)
					nj = f.size() - 1;
				else
					nj = j - 1;
				k2 = midpts[nj];
				/*quad face*/
				Facet fn;
				fn.push_back(f[j]);
				fn.push_back(k2);
				fn.push_back(fci);
				fn.push_back(k1);
				gFacets.push_back(fn);
			}
		/*triangles*/
		} else {
			/*add facets*/
			Int k1,k2;
			forEach(f,j) {
				k1 = f[j];
				nj = j + 1;
				if(nj == f.size())
					nj = 0;
				k2 = f[nj];
				/*quad face*/
				Facet fn;
				fn.push_back(k1);
				fn.push_back(fci);
				fn.push_back(k2);
				gFacets.push_back(fn);
			}
		}
	}
	/*************************
	 * Refine cells
	 *************************/
	forEach(rCells,i) {
		Int ci = rCells[i];
		Cell c = gCells[ci];
		Vector C = cC[ci];
		gVertices.push_back(C);
		Int cci = gVertices.size() - 1;
		Int ifBegin = gFacets.size();
		Cells newc;
		
		forEach(c,j) {
			Int fi = c[j];
			
			/*cell is refined but face is not?*/
			IntVector list;
			if(!refineF[fi]) {
				list.push_back(fi);
			} else {
				forEach(gFacets[fi],k)
					list.push_back(startF[fi] + k);
			}

			/*refine cells*/
			forEach(list,k) {
				Int fni = list[k];
				Facet f = gFacets[fni];
			
				Cell cn;
				cn.push_back(fni);

				Int v1i,v2i,nj;
				forEach(f,l) {
					v1i = f[l];
					nj = l + 1;
					if(nj == f.size())
						nj = 0;
					v2i = f[nj];
					//triangular face
					Facet fn;
					fn.push_back(v1i);
					fn.push_back(v2i);
					fn.push_back(cci);
					//duplicate
					Int k = ifBegin;
					for(; k < gFacets.size();k++) {
						if(equal(fn,gFacets[k])) {
							cn.push_back(k);
							break;
						}
					}
					if(k == gFacets.size()) {
						gFacets.push_back(fn);
						cn.push_back(k);
					}
				}
				newc.push_back(cn);
			}
		}
		/*merge cells*/
		if(rparams.shape == 0) {
			/*cells on same face*/
			IntVector ownerf;
			forEach(c,j) {
				Int fi = c[j];
				forEach(gFacets[fi],k)
					ownerf.push_back(fi);
			}
			/*find cells to merge with a shared face*/
			Cells mergec;
			forEach(newc,j) {
				Int o1 = ownerf[j];
				Cell& c1 = newc[j];
				IntVector mg;
				forEachS(newc,k,j+1) {
					Int o2 = ownerf[k];
					Cell& c2 = newc[k];
					if(o1 == o2) continue;
					forEach(c1,m) {
						Int f1 = c1[m];
						forEach(c2,n) {
							Int f2 = c2[n];
							if(f1 == f2) {
								mg.push_back(k);
								c1.erase(c1.begin() + m);
								c2.erase(c2.begin() + n);
								rsFacets.push_back(f1);
								goto END;
							}
						}
					}
					END:;					
				}
				mergec.push_back(mg);
			}
			/*merge cells*/
			IntVector erasei;
			erasei.assign(newc.size(),0);
			forEachR(mergec,j) {
				IntVector& cm = mergec[j];
				forEach(cm,k)
					erasei[cm[k]] = 1;
				if(cm.size()) {
					Cell& c1 = newc[j];
					Cell& c2 = newc[cm[0]];
					c1.insert(c1.end(),c2.begin(),c2.end());
				}
			}
			/*erase cells*/
			IntVector erasec;
			forEach(erasei,j) {
				if(erasei[j])
					erasec.push_back(j);
			}
			erase_indices(newc,erasec);
			/*two cells should share only one face*/
			forEach(newc,j) {
				Cell& c1 = newc[j];
				forEachS(newc,k,j+1) {
					Cell& c2 = newc[k];
					//find shared faces
					IntVector shared1,shared2;
					forEach(c1,m) {
						Int f1 = c1[m];
						forEach(c2,n) {
							Int f2 = c2[n];
							if(f1 == f2) {
								shared1.push_back(m);
								shared2.push_back(n);
							}
						}
					}
					//two faces shared between two cells
					if(shared1.size() == 2) {
						Int fi1 = c1[shared1[0]];
						Int fi2 = c1[shared1[1]];
						Facet& f1 = gFacets[fi1];
						Facet& f2 = gFacets[fi2];
						erase_indices(c1,shared1);
						erase_indices(c2,shared2);
						rsFacets.push_back(fi1);
						rsFacets.push_back(fi2);
						//add new face by merging two faces
						Facet f;
						if(f1[0] == f2[0]) {
							f.push_back(f1[1]);
							f.push_back(f1[0]);
							f.push_back(f2[1]);
							f.push_back(f2[2]);
						} else if(f1[0] == f2[1]) {
							f.push_back(f1[1]);
							f.push_back(f1[0]);
							f.push_back(f2[0]);
							f.push_back(f2[2]);
						} else if(f1[1] == f2[1]) {
							f.push_back(f1[0]);
							f.push_back(f1[1]);
							f.push_back(f2[0]);
							f.push_back(f2[2]);
						} else if(f1[1] == f2[0]) {
							f.push_back(f1[0]);
							f.push_back(f1[1]);
							f.push_back(f2[1]);
							f.push_back(f2[2]);
						}
						gFacets.push_back(f);
						Int fi = gFacets.size() - 1;
						c1.push_back(fi);
						c2.push_back(fi);
					}
				}
			}
		}
		/*add cells*/
		gCells.insert(gCells.end(),newc.begin(),newc.end());
		forEach(newc,j)
			cellMap.push_back(ci);
	}
	/*************************************
	 *  Remove refined facets and cells
	 *************************************/
	/*add additional facets to remove*/
	refineF.resize(gFacets.size(),0);
	forEach(rsFacets,i) {
		Int fi = rsFacets[i];
		if(!refineF[fi]) {
			refineF[fi] = Constants::MAX_INT - 1;
			rFacets.push_back(fi);
		}
	}
	/*adjust facet indexes in cells and boundaries*/
	Int count = 0;
	forEach(refineF,i) {
		if(refineF[i] == 0) refineF[i] = count++;
		else if(refineF[i] == Constants::MAX_INT - 1);
		else refineF[i] = Constants::MAX_INT;
	}
	/*erase old facets and add the new ones*/
#define ERASEADD() {									\
	IntVector newF,eraseF;								\
	forEach(c,j) {										\
		Int fi = c[j];									\
		if(refineF[fi] == Constants::MAX_INT) {			\
			eraseF.push_back(j);						\
			forEach(gFacets[fi],k) {					\
				Int fni = startF[fi] + k;				\
				fni = refineF[fni];						\
				newF.push_back(fni);					\
			}											\
		}												\
		c[j] = refineF[fi];								\
	}													\
	erase_indices(c,eraseF);							\
	c.insert(c.end(),newF.begin(),newF.end());			\
}
	forEach(gCells,i) {
		Cell& c = gCells[i];
		ERASEADD();
	}
	forEachIt(Boundaries,gBoundaries,it) {
		IntVector& c = it->second;
		ERASEADD();
	}
#undef ERASEADD

	/*erase facets*/
	sort(rFacets.begin(),rFacets.end());
	erase_indices(gFacets,rFacets);
	
	/*erase cells*/
	sort(rCells.begin(),rCells.end());
	erase_indices(gCells,rCells);
	erase_indices(cellMap,rCells);
	gBCS = gCells.size();
	gBCSfield = gBCS * DG::NP;
	
	/*******************************************************
	 * Break edge of faces that are not set for refinement 
	 * but should be due to neighboring refined faces 
	 *******************************************************/
	if(rparams.shape == 0) {
		forEach(gFacets,i) {
			Facet& f = gFacets[i];
			IntVector addf;
			Vector v1,v2;
			addf.assign(f.size(),Constants::MAX_INT);
			forEach(f,j) {
				v1 = gVertices[f[j]];
				nj = j + 1;
				if(nj == f.size())
					nj = 0;
				v2 = gVertices[f[nj]];
				Vector Ce = (v1 + v2) / 2.0;
		
				//duplicate
				Int k = ivBegin;
				for(; k < gVertices.size();k++) {
					if(equal(Ce,gVertices[k])) {
						addf[j] = k;
						break;
					}
				}
			}
			if(addf.size()) {
				Facet nf;
				forEach(f,j) {
					nf.push_back(f[j]);
					if(addf[j] != Constants::MAX_INT)
						nf.push_back(addf[j]);
				}
				f = nf;
			}
		}
	}
}
/**
create fields
*/
void Prepare::createFields(vector<string>& fields,Int step) {
	BaseField::destroyFields();

	/*for each field*/
	forEach(fields,i) {
		/*read at time 0*/
		stringstream path;
		Int size;
		path << fields[i] << step;
		std::string str = path.str(); 

		ifstream is(str.c_str());
		if(!is.fail()) {
			/*fields*/
			is >> str >> size;
			BaseField* bf;
			switch(size) {
				case 1 :  bf = new ScalarCellField(fields[i].c_str(),READWRITE,false); break;
				case 3 :  bf = new VectorCellField(fields[i].c_str(),READWRITE,false); break;
				case 6 :  bf = new STensorCellField(fields[i].c_str(),READWRITE,false); break;
				case 9 :  bf = new TensorCellField(fields[i].c_str(),READWRITE,false); break;
			}
			/*end*/
		}
	}
}
/**
open fields
*/
Int Prepare::readFields(vector<string>& fields,Int step) {
	Int count = 0;
	forEach(fields,i) {
		stringstream fpath;
		fpath << fields[i] << step;
		ifstream is(fpath.str().c_str());
		if(is.fail())
			continue;
		count++;
		break;
	}
	if(count)
		Mesh::read_fields(step);
	return count;
}
