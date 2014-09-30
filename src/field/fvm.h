#ifndef __FVM_H
#define __FVM_H

/* ***************************************
 * Implicit boundary conditions
 * ***************************************/
 template <class T> 
 void applyImplicitBCs(const MeshMatrix<T>& M) {
	 using namespace Mesh;
	 MeshField<T,CELL>& cF = *M.cF;
	 BasicBCondition* bbc;
	 BCondition<T>* bc;

	 /*boundary conditions*/
	 forEach(AllBConditions,i) {
		 bbc = AllBConditions[i];
		 if(bbc->fIndex == cF.fIndex) {
			 if(bbc->cIndex == NEUMANN ||
				 bbc->cIndex == SYMMETRY)
				 ;
			 else continue;

			 bc = static_cast<BCondition<T>*> (bbc);
			 Int sz = bc->bdry->size();
			 if(sz == 0) continue;

			 for(Int j = 0;j < sz;j++) {
				 Int k = (*bc->bdry)[j];
				 Int c1 = gFO[k];
				 Int c2 = gFN[k];
				 if(bc->cIndex == NEUMANN) {
					 Vector dv = cC[c2] - cC[c1];
					 M.ap[c1] -= M.an[1][k];
					 M.Su[c1] += M.an[1][k] * (bc->value * mag(dv));
					 M.an[1][k] = 0;
				 } else if(bc->cIndex == ROBIN) {
					 Vector dv = cC[c2] - cC[c1];
					 M.ap[c1] -= (1 - bc->shape) * M.an[1][k];
					 M.Su[c1] += M.an[1][k] * (bc->shape * bc->value + 
						 (1 - bc->shape) * bc->tvalue * mag(dv));
					 M.an[1][k] = 0;
				 } else if(bc->cIndex == SYMMETRY) {
					 M.ap[c1] -= M.an[1][k];
					 M.Su[c1] += M.an[1][k] * (sym(cF[c1],fN[k]) - cF[c1]);
					 M.an[1][k] = 0;
				 }
			 }
		 }
	 }
 }
 /*************************************
  * Exchange ghost cell information
  *************************************/
 template <class T> 
 void exchange_ghost(T* P) {
 	using namespace Mesh;
 	/*blocked exchange*/
 	if(Controls::ghost_exchange == Controls::BLOCKED) {
 		MeshField<T,CELL> buffer;
 		forEach(gInterMesh,i) {
 			interBoundary& b = gInterMesh[i];
 			IntVector& f = *(b.f);
 			if(b.from < b.to) {
 				//send
 				forEach(f,j)
 					buffer[j] = P[gFO[f[j]]];
 				MP::send(&buffer[0],f.size(),b.to,MP::FIELD);
 				//recieve
 				MP::recieve(&buffer[0],f.size(),b.to,MP::FIELD);
 				forEach(f,j)
 					P[gFN[f[j]]] = buffer[j];
 			} else {
 				//recieve
 				MP::recieve(&buffer[0],f.size(),b.to,MP::FIELD);
 				forEach(f,j)
 					P[gFN[f[j]]] = buffer[j];
 				//send 
 				forEach(f,j)
 					buffer[j] = P[gFO[f[j]]];
 				MP::send(&buffer[0],f.size(),b.to,MP::FIELD);
 			}
 		}
    /*Asynchronous exchange*/
 	} else {
 		MeshField<T,CELL> sendbuf,recvbuf;
 		std::vector<MP::REQUEST> request(2 * gInterMesh.size(),0);
 		Int rcount = 0;
 		//fill send buffer
 		forEach(gInterMesh,i) {
 			interBoundary& b = gInterMesh[i];
 			IntVector& f = *(b.f);
 			forEach(f,j) 
 				sendbuf[b.buffer_index + j] = P[gFO[f[j]]];
 		}

 		forEach(gInterMesh,i) {
 			interBoundary& b = gInterMesh[i];
 			//non-blocking send/recive
 			MP::isend(&sendbuf[b.buffer_index],b.f->size(),
 				b.to,MP::FIELD,&request[rcount]);
 			rcount++;
 			MP::irecieve(&recvbuf[b.buffer_index],b.f->size(),
 				b.to,MP::FIELD,&request[rcount]);
 			rcount++;
 		}
 		//wait
 		MP::waitall(rcount,&request[0]);
 		//recieve buffer
 		forEach(gInterMesh,i) {
 			interBoundary& b = gInterMesh[i];
 			IntVector& f = *(b.f);
 			forEach(f,j)
 				P[gFN[f[j]]] = recvbuf[b.buffer_index + j];
 		}
 	}
 	/*end*/
 }

/* ***************************************
 * Explicit boundary conditions
 * **************************************/
template<class T,ENTITY E>
void updateExplicitBCs(const MeshField<T,E>& cF,
							  bool update_ghost = false,
							  bool update_fixed = false
							  ) {
	using namespace Mesh;
	BasicBCondition* bbc;
	BCondition<T>* bc;
	Scalar z = Scalar(0),
		   zmin = Scalar(0),
		   zmax = Scalar(0),
		   zR = Scalar(0);
	Vector C(0);

	/*boundary conditions*/
	forEach(AllBConditions,i) {
		bbc = AllBConditions[i];
		if(bbc->fIndex == cF.fIndex) {
			if(bbc->cIndex == GHOST) 
				continue;

			bc = static_cast<BCondition<T>*> (bbc);
			Int sz = bc->bdry->size();
			if(sz == 0) continue;

			if(update_fixed) {
				if(bc->cIndex == DIRICHLET || 
					bc->cIndex == POWER || 
					bc->cIndex == LOG || 
					bc->cIndex == PARABOLIC ||
					bc->cIndex == INVERSE
					) {
						Int ci,j;
						Scalar r;
						if(bc->zMax > 0) {
							zmin = bc->zMin;
							zmax = bc->zMax;
							zR = zmax - zmin;
						} else {
							zmin = Scalar(10e30);
							zmax = -Scalar(10e30);
							C = Vector(0);
							for(j = 0;j < sz;j++) {
								Facet& f = gFacets[j];
								forEach(f,k) {
									z = (vC[f[k]] & bc->dir);
									if(z < zmin) 
										zmin = z;
									if(z > zmax) 
										zmax = z;
								}
								C += fC[j];
							}
							C /= Scalar(sz);
							zR = zmax - zmin;

							if(bc->cIndex == PARABOLIC) {
								ci = gFN[(*bc->bdry)[0]];
								zR = magSq(cC[ci] - C);
								for(j = 1;j < sz;j++) {
									ci = gFN[(*bc->bdry)[0]];
									r = magSq(cC[ci] - C);
									if(r < zR) zR = r;
								}
							}
						}
				}
			}
			for(Int j = 0;j < sz;j++) {
				Int k = (*bc->bdry)[j];
				Int c1 = gFO[k];
				Int c2 = gFN[k];
				if(bc->cIndex == NEUMANN) {
					Vector dv = cC[c2] - cC[c1];
					cF[c2] = cF[c1] + bc->value * mag(dv);
				} else if(bc->cIndex == ROBIN) {
					Vector dv = cC[c2] - cC[c1];
					cF[c2] = bc->shape * bc->value + 
						(1 - bc->shape) * (cF[c1] + bc->tvalue * mag(dv));
				} else if(bc->cIndex == SYMMETRY) {
					cF[c2] = sym(cF[c1],fN[k]);
				} else if(bc->cIndex == CYCLIC) {
					Int c22;
					if(j < sz / 2) 
						c22 = gFO[(*bc->bdry)[j + sz/2]];
					else
						c22 = gFO[(*bc->bdry)[j - sz/2]];
					cF[c2] = cF[c22];
				} else { 
					if(update_fixed) {
						T v(0);
						z = (cC[c2] & bc->dir) - zmin;
						if(bc->cIndex == DIRICHLET) {
							v = bc->value;
						} else if(bc->cIndex == POWER) {
							if(z < 0) z = 0;
							if(z > zR) v = bc->value;
							else v = bc->value * pow(z / zR,bc->shape);
						} else if(bc->cIndex == LOG) {
							if(z < 0) z = 0;
							if(z > zR) v = bc->value;
							else v = bc->value * (log(1 + z / bc->shape) / log(1 + zR / bc->shape));
						} else if(bc->cIndex == PARABOLIC) {
							z = magSq(cC[c2] - C);
							v = bc->value * (z / zR);
						} else if(bc->cIndex == INVERSE) {
							v = bc->value / (z + bc->shape);
						}
						if(!bc->first && !equal(mag(bc->tvalue),0)) { 
							T meanTI = v * (bc->tvalue * pow (z / zR,-bc->tshape));
							Scalar rFactor = 4 * ((rand() / Scalar(RAND_MAX)) - 0.5);
							v += ((cF[c2] - v) * 0.9 + (meanTI * rFactor) * 0.1);
						}
						bc->fixed[j] = cF[c2] = v;
					} else {
						cF[c2] = bc->fixed[j];
					}
				}
			}
			bc->first = false;
		}
	}
	/*ghost cells*/
	if(update_ghost && gInterMesh.size()) {
		exchange_ghost(&cF[0]);
	}
}
/* ***************************************
 * Fill boundary from internal values
 * **************************************/
template<class T,ENTITY E>
const MeshField<T,E>& fillBCs(const MeshField<T,E>& cF,
			 bool update_ghost = false) {
	/*neumann update*/
	using namespace Mesh;
	forEachS(cF,i,gBCellsStart)
		cF[i] = cF[gFO[gCells[i][0]]];
	/*ghost cells*/
	if(update_ghost && gInterMesh.size()) {
		exchange_ghost(&cF[0]);
	}
	return cF;
}
/* *******************************
 * matrix - vector product p * q
 * *******************************/
template <class T> 
MeshField<T,CELL> mul (const MeshMatrix<T>& p,const MeshField<T,CELL>& q) {
	using namespace Mesh;
	MeshField<T,CELL> r;
	Int c1,c2;
	r = q * p.ap;
	forEach(gFacets,f) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] -= q[c2] * p.an[1][f];
		r[c2] -= q[c1] * p.an[0][f];
	}
	return r;
}
/*matrix transopose - vector product pT * q */
template <class T> 
MeshField<T,CELL> mult (const MeshMatrix<T>& p,const MeshField<T,CELL>& q) {
	using namespace Mesh;
	MeshField<T,CELL> r;
	Int c1,c2;
	r = q * p.ap;
	forEach(gFacets,f) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c2] -= q[c1] * p.an[1][f];
		r[c1] -= q[c2] * p.an[0][f];
	}
	return r;
}
/* calculate RHS sum */
template <class T> 
MeshField<T,CELL> getRHS(const MeshMatrix<T>& p) {
	using namespace Mesh;
	MeshField<T,CELL> r;
	Int c1,c2;
	r = p.Su;
	forEach(gFacets,f) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] += (*p.cF)[c2] * p.an[1][f];
		r[c2] += (*p.cF)[c1] * p.an[0][f];
	}
	return r;
}

/* ********************************
 * Interpolate field operations
 * *******************************/
/*central difference*/
template<class type>
MeshField<type,FACET> cds(const MeshField<type,CELL>& cF) {
	using namespace Mesh;
	MeshField<type,FACET> fF;
	forEach(fF,i) {
		fF[i] =  (cF[gFO[i]] * (fI[i])) + (cF[gFN[i]] * (1 - fI[i]));
	}
	return fF;
}
/*upwind*/
template<class type>
MeshField<type,FACET> uds(const MeshField<type,CELL>& cF,const ScalarFacetField& flux) {
	using namespace Mesh;
	MeshField<type,FACET> fF;
	forEach(fF,i) {
		if(flux[i] >= 0) fF[i] = cF[gFO[i]];
		else fF[i] = cF[gFN[i]];
	}
	return fF;
}
/*facet data to vertex data */
template<class type>
MeshField<type,VERTEX> cds(const MeshField<type,FACET>& fF) {
	using namespace Mesh;
	std::vector<Scalar> cnt;
	MeshField<type,VERTEX> vF;
	cnt.assign(vF.size(),Scalar(0));
	Scalar dist;

	vF = type(0);
	forEach(fF,i) {
		Facet& f = gFacets[i];
		if(gFN[i] < gBCellsStart) {
			forEach(f,j) {
				dist = 1.f / magSq(gVertices[f[j]] - fC[i]);
				vF[f[j]] += (fF[i] * dist);
				cnt[f[j]] += dist;
			}
		} else {
			forEach(f,j) {
				vF[f[j]] += Scalar(10e30) * fF[i];
				cnt[f[j]] += Scalar(10e30);
			}
		}
	}
	forEach(vF,i) {
		vF[i] /= cnt[i];
		if(mag(vF[i]) < Constants::MachineEpsilon) 
			vF[i] = type(0);
	}
	return vF;
}
/* **************************
 * Integrate field operation
 * **************************/
template<class type>
MeshField<type,CELL> sum(const MeshField<type,FACET>& fF) {
	using namespace Mesh;
	MeshField<type,CELL> cF;
	cF = type(0);
	forEach(fF,i) {
		cF[gFO[i]] += fF[i];
		cF[gFN[i]] -= fF[i];
	}
	return cF;
}
/*****************************************
 * Non-integrated field operations. The field operations to be defined later
 * such as div,lap,grad,ddt,ddt2,src etc.. are values integrated over a volume. 
 * The following macro versions give the corresponding non-integrated cell
 * center values by dividing with the cell volumes
 *****************************************/
#define divi(x)  fillBCs((div(x)   / Mesh::cV),true)
#define lapi(x)  fillBCs((lap(x)   / Mesh::cV),true)
#define ddti(x)  fillBCs((ddt(x)   / Mesh::cV),true)
#define ddt2i(x) fillBCs((ddt2(x)  / Mesh::cV),true)
#define srci(x)  fillBCs((src(x)   / Mesh::cV),true)
#define gradi(x) fillBCs((grad(x)  / Mesh::cV),true)

/**********************************************************************
 * Gradient field operation.
 **********************************************************************/

/*Explicit*/
inline VectorCellField grad(const ScalarFacetField& p) {
	return sum(mul(Mesh::fN,p));
}
inline VectorCellField grad(const ScalarCellField& p) {
	return grad(cds(p));
}
inline TensorCellField grad(const VectorFacetField& p) {
	return sum(mul(Mesh::fN,p));
}
inline TensorCellField grad(const VectorCellField& p) {
	return grad(cds(p));
}

/* *********************************************
 * Laplacian field operation
 * ********************************************/

/*Implicit*/
template<class type>
MeshMatrix<type> lap(MeshField<type,CELL>& cF,const ScalarFacetField& mu) {
	using namespace Controls;
	using namespace Mesh;
	MeshMatrix<type> m;
	VectorFacetField K;
	Vector dv;
	Int c1,c2;
	Scalar D = 0;
	/*clear*/
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	m.Su = type(0);
	m.ap = Scalar(0);
	forEach(mu,i) {
		c1 = gFO[i];
		c2 = gFN[i];
		dv = cC[c2] - cC[c1];
		/*diffusivity coefficient*/
		if(nonortho_scheme == NONE) {
			D = mag(fN[i]) / mag(dv);
		} else {
			if(nonortho_scheme == OVER_RELAXED) {
				D = ((fN[i] & fN[i]) / (fN[i] & dv));
			} else if(nonortho_scheme == MINIMUM) {
				D = ((fN[i] & dv) / (dv & dv));
			} else if(nonortho_scheme == ORTHOGONAL) {
				D = sqrt((fN[i] & fN[i]) / (dv & dv));
			}
			K[i] = fN[i] - D * dv;
		}
		/*coefficients*/
		m.an[0][i] = D * mu[i];
		m.an[1][i] = D * mu[i];
		m.ap[c1]  += m.an[0][i];
		m.ap[c2]  += m.an[1][i];
	}
	/*non-orthogonality handled through deferred correction*/
	if(nonortho_scheme != NONE) {
		MeshField<type,FACET> r = dot(cds(gradi(cF)),K);
		type res;
		forEach(mu,i) {
			c1 = gFO[i];
			c2 = gFN[i];
			res = m.an[0][i] * (cF[c2] - cF[c1]);
			if(mag(r[i]) > Scalar(0.5) * mag(res)) 
				r[i] = Scalar(0.5) * res;
		}
		m.Su = sum(r);
	}
	/*end*/
	return m;
}

template<class type>
inline MeshMatrix<type> lap(MeshField<type,CELL>& cF,const ScalarCellField& mu) {
	return lap(cF,cds(mu));
}

/* ***************************************************
 * Divergence field operation
 * ***************************************************/ 
/*face flux*/
inline ScalarFacetField flx(const VectorFacetField& p) {
	return dot(p,Mesh::fN);
}
inline ScalarFacetField flx(const VectorCellField& p) {
	return flx(cds(p));
}
inline VectorFacetField flx(const TensorFacetField& p) {
	return dot(p,Mesh::fN);
}
inline VectorFacetField flx(const TensorCellField& p) {
	return flx(cds(p));
}
/* Explicit */
inline ScalarCellField div(const VectorFacetField& p) {
	return sum(flx(p));
}
inline ScalarCellField div(const VectorCellField& p) {
	return sum(flx(p));
}
inline VectorCellField div(const TensorFacetField& p) {
	return sum(flx(p));
}
inline VectorCellField div(const TensorCellField& p) {
	return sum(flx(p));
}
/* Implicit */
template<class type>
MeshMatrix<type> div(MeshField<type,CELL>& cF,const ScalarFacetField& flux,const ScalarFacetField& mu) {
	using namespace Controls;
	using namespace Mesh;
	MeshMatrix<type> m;
	Scalar F,G;
	m.cF = &cF;
	m.flags = 0;
	m.Su = type(0);
	m.ap = Scalar(0);

	/*Implicit convection schemes*/
	bool isImplicit = (
		convection_scheme == CDS ||
		convection_scheme == UDS ||
		convection_scheme == BLENDED ||
		convection_scheme == HYBRID );

	if(isImplicit) {
		ScalarFacetField gamma;
		if(convection_scheme == CDS) 
			gamma = Scalar(1);
		else if(convection_scheme == UDS) 
			gamma = Scalar(0);
		else if(convection_scheme == BLENDED) 
			gamma = Scalar(blend_factor);
		else if(convection_scheme == HYBRID) {
			Scalar D;
			Vector dv;
			forEach(gFacets,j) {
				/*calc D - uncorrected */
				dv = cC[gFN[j]] - cC[gFO[j]];
				D = (mag(fN[j]) / mag(dv)) * mu[j];
				/*compare F and D */
				F = flux[j];
				if(F < 0) {
					if(-F * fI[j] > D) gamma[j] = 0;
					else gamma[j] = 1;
				} else {
					if(F * (1 - fI[j]) > D) gamma[j] = 0;
					else gamma[j] = 1;
				}
			}
		}
		forEach(flux,i) {
			F = flux[i];
			G = gamma[i];
			m.an[0][i] = ((G) * (-F * (  fI[i]  )) + (1 - G) * (-max( F,0)));
			m.an[1][i] = ((G) * ( F * (1 - fI[i])) + (1 - G) * (-max(-F,0)));
			m.ap[gFO[i]] += m.an[0][i];
			m.ap[gFN[i]] += m.an[1][i];
		}
	/*deferred correction*/
	} else {
		forEach(flux,i) {
			F = flux[i];
			m.an[0][i] = -max( F,0);
			m.an[1][i] = -max(-F,0);
			m.ap[gFO[i]] += m.an[0][i];
			m.ap[gFN[i]] += m.an[1][i];
		}

		MeshField<type,FACET> corr;
		if(convection_scheme == CDSS) {
			corr = cds(cF) - uds(cF,flux);
		} else if(convection_scheme == LUD) {
			VectorFacetField R = fC - uds(cC,flux);
			corr = dot(uds(gradi(cF),flux),R);
		} else if(convection_scheme == MUSCL) {
			VectorFacetField R = fC - uds(cC,flux);
			corr  = (  blend_factor  ) * (cds(cF) - uds(cF,flux));
			corr += (1 - blend_factor) * (dot(uds(gradi(cF),flux),R));
		} else {
			/*
			TVD schemes
			~~~~~~~~~~~
			Reference:
				M.S Darwish and F Moukalled "TVD schemes for unstructured grids"
				Versteeg and Malaskara
			Description:
				phi = phiU + psi(r) * [(phiD - phiC) * (1 - fi)]
			Schemes
				psi(r) = 0 =>UDS
				psi(r) = 1 =>CDS
			R is calculated as ratio of upwind and downwind gradient
				r = phiDC / phiCU
			Further modification to unstructured grid to better fit LUD scheme
				r = (phiDC / phiCU) * (fi / (1 - fi))
			*/
			/*calculate r*/
			MeshField<type,FACET> q,r,phiDC,phiCU;
			ScalarFacetField uFI;
			{
				ScalarFacetField nflux = Scalar(0)-flux;
				phiDC = uds(cF,nflux) - uds(cF,flux);
				forEach(phiDC,i) {
					if(flux[i] >= 0) G = fI[i];
					else G = 1 - fI[i];
					uFI[i] = G;
				}
				/*Bruner's or Darwish way of calculating r*/
				if(TVDbruner) {
					VectorFacetField R = fC - uds(cC,flux);
					phiCU = 2 * (dot(uds(gradi(cF),flux),R));
				} else {
					VectorFacetField R = uds(cC,nflux) - uds(cC,flux);
					phiCU = 2 * (dot(uds(gradi(cF),flux),R)) - phiDC;
				}
				/*end*/
			}
			r = (phiCU / phiDC) * (uFI / (1 - uFI));
			forEach(phiDC,i) {
				if(equal(phiDC[i] * (1 - uFI[i]),type(0)))
					r[i] = type(0);
			}
			/*TVD schemes*/
			if(convection_scheme == VANLEER) {
				q = (r+fabs(r)) / (1+r);
			} else if(convection_scheme == VANALBADA) {
				q = (r+r*r) / (1+r*r);
			} else if(convection_scheme == MINMOD) {
				q = max(type(0),min(r,type(1)));
			} else if(convection_scheme == SUPERBEE) {
				q = max(min(r,type(2)),min(2*r,type(1)));
				q = max(q,type(0));
			} else if(convection_scheme == SWEBY) {
				Scalar beta = 2;
				q = max(min(r,type(beta)),min(beta*r,type(1)));
				q = max(q,type(0));
			} else if(convection_scheme == QUICKL) {
				q = min(2*r,(3+r)/4);
				q = min(q,type(2));
				q = max(q,type(0));
			} else if(convection_scheme == UMIST) {
				q = min(2*r,(3+r)/4);
				q = min(q,(1+3*r)/4);
				q = min(q,type(2));
				q = max(q,type(0));
			} else if(convection_scheme == QUICK) {
				q = (3+r)/4;
			} else if(convection_scheme == DDS) {
				q = 2;
			} else if(convection_scheme == FROMM) {
				q = (1+r)/2;
			}
			corr = q * phiDC * (1 - uFI);
			/*end*/
		}
		m.Su = sum(flux * corr);
	}
	return m;
}

template<class type,ENTITY E>
inline MeshMatrix<type> div(MeshField<type,CELL>& cF,const MeshField<Vector,E>& rhoU,const ScalarFacetField& mu) {
	return div(cF,div(rhoU),mu);
}

/* *******************************
 * Temporal derivative
 * *******************************/
template<class type>
MeshMatrix<type> ddt(MeshField<type,CELL>& cF,const ScalarCellField& rho) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	if(Controls::time_scheme == Controls::EULER || !(cF.access & STOREPREV)) {
		if(Controls::time_scheme != Controls::EULER) cF.initStore();
		m.ap = (Mesh::cV * rho) / -Controls::dt;
		m.Su = cF * m.ap;
	} else if(Controls::time_scheme == Controls::SECOND_ORDER) {
		m.ap = (1.5 * Mesh::cV * rho) / -Controls::dt;
		m.Su = ((4.0 * cF - cF.tstore[1]) / 3.0) * m.ap;
	}
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	return m;
}
template<class type>
MeshMatrix<type> ddt2(MeshField<type,CELL>& cF,const ScalarCellField& rho) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	if(!(cF.access & STOREPREV)) cF.initStore();
	m.ap = (Mesh::cV * rho) / -(Controls::dt * Controls::dt);
	m.Su = (2.0 * cF - cF.tstore[1]) * m.ap;
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	return m;
}
/* *******************************
 * Linearized source term
 * *******************************/
template<class type>
MeshMatrix<type> src(MeshField<type,CELL>& cF,const MeshField<type,CELL>& Su,const ScalarCellField& Sp) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	m.ap = -(Sp * Mesh::cV);
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	m.Su = (Su * Mesh::cV);
	return m;
}
/*Explicit*/
template<class type>
MeshField<type,CELL> src(const MeshField<type,CELL>& Su) {
	return (Su * Mesh::cV);
}
/* ************************************************
 * Form transport equation
 *************************************************/
template<class type>
void addTemporal(MeshMatrix<type>& M,const ScalarCellField& rho,Scalar cF_UR) {
	using namespace Controls;
	if(state == STEADY)
		M.Relax(cF_UR);
	else {
		if(!equal(implicit_factor,1)) {
			MeshField<type,CELL> k1 = M.Su - mul(M, *M.cF);
			if(runge_kutta == 1) {
				M = M * (implicit_factor) + k1  * (1 - implicit_factor);
			} else {
				ScalarCellField mdt = Controls::dt / (rho * Mesh::cV);
				if (runge_kutta == 2) {
					MeshField<type, CELL> k2 = M.Su - mul(M, *M.cF + k1 * mdt);
					M = (M * (implicit_factor) + k2 * (1 - implicit_factor) + k1) / 2;
				} else if (runge_kutta == 3) {
					MeshField<type, CELL> k2 = M.Su - mul(M, *M.cF + k1 * mdt / 2);
					MeshField<type, CELL> k3 = M.Su - mul(M, *M.cF + (2 * k2 - k1) * mdt);
					M = (M * (implicit_factor) + k3 * (1 - implicit_factor) + 4 * k2 + k1) / 6;
				} else if (runge_kutta == 4) {
					MeshField<type, CELL> k2 = M.Su - mul(M, *M.cF + k1 * mdt / 2);
					MeshField<type, CELL> k3 = M.Su - mul(M, *M.cF + k2 * mdt / 2);
					MeshField<type, CELL> k4 = M.Su - mul(M, *M.cF + k3 * mdt);
					M = (M * (implicit_factor) + k4 * (1 - implicit_factor) + 2 * k3 + 2 * k2 + k1) / 6;
				}
			}
			if(equal(implicit_factor,0))
				M.flags = M.SYMMETRIC;
		}
		M += ddt(*M.cF,rho);
	}
}

template<class type>
MeshMatrix<type> diffusion(MeshField<type,CELL>& cF,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR) {
	MeshMatrix<type> M = -lap(cF,mu);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> convection(MeshField<type,CELL>& cF,const ScalarFacetField& F,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR) {
	MeshMatrix<type> M = div(cF,F,mu);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> transport(MeshField<type,CELL>& cF,const ScalarFacetField& F,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR) {
	MeshMatrix<type> M = div(cF,F,mu) - lap(cF,mu);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> transport(MeshField<type,CELL>& cF,const ScalarFacetField& F,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR, 
		const MeshField<type,CELL>& Su,const ScalarCellField& Sp) {
	MeshMatrix<type> M = div(cF,F,mu) - lap(cF,mu) - src(cF,Su,Sp);
	addTemporal(M,rho,cF_UR);
	return M;
}

#endif