#include "mesh.h"

using namespace std;

/*global mesh*/
namespace Mesh {
	MeshObject        gMesh;
	std::string&      gMeshName = gMesh.name;
	Vertices&         gVertices = gMesh.mVertices;
	Facets&           gFacets   = gMesh.mFacets;
	Cells&            gCells    = gMesh.mCells;
	Boundaries&       gBoundaries = gMesh.mBoundaries;
	IntVector&        gFOC = gMesh.mFOC;
	IntVector&        gFNC = gMesh.mFNC;
	Int&              gBCS = gMesh.mBCS;
	Int&              gBCSI = gMesh.mBCSI;
	Cells&            gFaceID = gMesh.mFaceID;
	InterBoundVector& gInterMesh = gMesh.mInterMesh;
	VectorVector&     gfC = gMesh.mFC;
	VectorVector&     gcC = gMesh.mCC;
	VectorVector&     gfN = gMesh.mFN;
	ScalarVector&     gcV = gMesh.mCV;
	
	vector<BasicBCondition*> AllBConditions;
	Vertices         probePoints;
}
namespace Controls {
	RefineParams refine_params;
}
void Controls::enrollRefine(Util::ParamList& params) {
	params.enroll("direction",&refine_params.dir);
	params.enroll("field",&refine_params.field);
	params.enroll("field_max",&refine_params.field_max);
	params.enroll("field_min",&refine_params.field_min);
	params.enroll("limit",&refine_params.limit);
}
void Mesh::clear() {
	gMesh.clear();
	Mesh::clearBC();
	probePoints.clear();
}
void Mesh::MeshObject::clear() {
	mVertices.clear();
	mFacets.clear();
	mCells.clear();
	mBoundaries.clear();
	mFOC.clear();
	mFNC.clear();
	mInterMesh.clear();
	mPatches.clear();
}
/*read mesh*/
bool Mesh::MeshObject::readMesh(Int step,bool first, bool coarse) {
	/*open file*/
	stringstream path;
	if(coarse) 
		path << name;
	else
		path << name << "_" << step;
	string str = path.str();
	ifstream is(str.c_str());
	if(is.fail()) {
		if(first) {
			str = name;
			is.open(str.c_str());
		} else
			return false;
	}
	/*read*/
	clear();
	is >> hex;
	is >> mVertices;
	is >> mFacets;
	is >> mCells;
	while(Util::nextc(is)) {
		IntVector index;
		string str;
		is >> str;
		is >> index;

		IntVector& gB = mBoundaries[str];
		gB.insert(gB.begin(),index.begin(),index.end());

		/*internal mesh boundaries*/
		if(str.find("interMesh") != std::string::npos) {
			interBoundary b;
			sscanf(str.c_str(), "interMesh_%x_%x", &b.from,&b.to);
			b.f    = &mBoundaries[str];
			mInterMesh.push_back(b);
		}
	}
	/*start of buffer*/ 
	Int buffer_index = 0;
	forEach(mInterMesh,i) {
		interBoundary& b = mInterMesh[i];
		b.buffer_index = buffer_index;
		buffer_index += b.f->size();
	}
	is >> dec;
	return true;
}
/*write mesh*/
void Mesh::MeshObject::writeMesh(ostream& os) {
	os << hex;
	os.precision(12);
	os << mVertices;
	os.precision(6);
	os << mFacets;
	os << mCells;
	forEachIt(Boundaries,mBoundaries,it)
		os << it->first << " " << it->second << endl;
	os << dec;
}
/*add boundary cells*/
void Mesh::MeshObject::addBoundaryCells() {
	using namespace Constants;
	
	/*neighbor and owner cells of face*/
	mBCS = mCells.size();
	mFOC.assign(mFacets.size(),MAX_INT);
	mFNC.assign(mFacets.size(),MAX_INT);
	forEach(mCells,i) {
		Cell& c = mCells[i];
		forEach(c,j) {
			Int fi = c[j];
			if(mFOC[fi] == MAX_INT) 
				mFOC[fi] = i;
			else 
				mFNC[fi] = i;
		}
	}
	/*Flag boundary faces not in mBoundaries for auto deletion*/
	{
		IntVector faceInB;
		faceInB.assign(mFacets.size(),0);
		forEachIt(Boundaries,mBoundaries,it) {
			IntVector& mB = it->second;	
			forEach(mB,j)
				faceInB[mB[j]] = 1;
		}
	
		IntVector& mDelete = mBoundaries["delete"];
		forEach(mFNC,i) {
			if(mFNC[i] == MAX_INT) {
				if(!faceInB[i])
					mDelete.push_back(i);
			}
		}
	}
	/*reorder cells*/
	{	
		Cells bcs;
		IntVector allbs;
		Int count = 0;
		Int bdry_size = 0;
		
		forEachIt(Boundaries,mBoundaries,it)
			bdry_size += it->second.size();
		bcs.resize(bdry_size);
		allbs.resize(bdry_size);
		forEach(mCells,i) {
			Cell& c = mCells[i];
			forEach(c,j) {
				Int fi = c[j];
				if(mFNC[fi] == MAX_INT) {
					allbs[count] = i;
					bcs[count] = c;
					count++;
					break;
				}
			}
		}
		mBCSI = mBCS - count;
		allbs.resize(count);
		bcs.resize(count);
		erase_indices(mCells,allbs);
		mCells.insert(mCells.end(),bcs.begin(),bcs.end());
	}
	/*add boundary cells*/
	forEachIt(Boundaries,mBoundaries,it) {
		IntVector& mB = it->second;
		forEach(mB,j) {
			Int fi = mB[j];
			if(mFNC[fi] == MAX_INT) {
				Cell c;
				c.push_back(fi);
				mCells.push_back(c);
				mFNC[fi] = mCells.size() - 1;
			}
		}
	}
}
void Mesh::MeshObject::calcGeometry() {
	Int i;
	/*allocate*/
	mFC.assign(mFacets.size(),Vector(0));
	mCC.assign(mCells.size(),Vector(0));
	mFN.assign(mFacets.size(),Vector(0));
	mCV.assign(mCells.size(),Scalar(0));
	mReversed.assign(mFacets.size(),false);

	/* face centre*/
	forEach(mFacets,i) {
		Facet& f = mFacets[i];
		Vector C(0);
		forEach(f,j)
			C += mVertices[f[j]];
		mFC[i] = C / Scalar(f.size());
	}
	/* cell centre */
	forEach(mCells,i) {
		Cell& c = mCells[i];
		Vector C(0);
		forEach(c,j)
			C += mFC[c[j]];
		mCC[i] = C / Scalar(c.size());
	}
	/* face normal */
	Vector v1,v2,v3,v;
	Scalar magN;
	forEach(mFacets,i) {
		Facet& f = mFacets[i];
		Vector N(0),C(0),Ci,Ni;
		Scalar Ntot = Scalar(0);
		v1 = mFC[i];
		forEach(f,j) {
			v2 = mVertices[f[j]];
			if(j + 1 == f.size())
				v3 = mVertices[f[0]];
			else
				v3 = mVertices[f[j + 1]];
			Ni = ((v2 - v1) ^ (v3 - v1));
			magN = mag(Ni);
			Ci = magN * ((v1 + v2 + v3) / 3);
			
			
			C += Ci;
			Ntot += magN;
			N += Ni;
		}
		mFC[i] = C / Ntot;    /*corrected face centre*/
		v = mFC[i] - mCC[mFOC[i]];
		if((v & N) < 0) {
			N = -N;
			mReversed[i] = true;
		}
		mFN[i] = N / Scalar(2);
	}
	/* cell volumes */
	for(i = 0;i < mBCS;i++) {
		Cell& c = mCells[i];
		Scalar V(0),Vi;
		Vector v = mCC[i],C(0);
		forEach(c,j) {
			v = mCC[i] - mFC[c[j]];
			Vi = mag(v & mFN[c[j]]);
			C += Vi * (3 * mFC[c[j]] + mCC[i]) / 4;
			V += Vi;
		}
		mCC[i] = C / V;         /*corrected cell centre */
		mCV[i] = V / Scalar(3);
	}
	/*boundary cell centre and volume*/
	forEachS(mCells,i,mBCS) {
		Int fi = mCells[i][0];
		mCV[i] = mCV[mFOC[fi]];
		mCC[i] = mFC[fi];
	}
	/*facet ids*/
	mFaceID.clear();
	forEach(mCells,i) {
		IntVector b;
		forEach(mCells[i],j)
			b.push_back(j);
		mFaceID.push_back(b);
	}
}
/* 
* Remove empty boundary
*/
void Mesh::MeshObject::removeBoundary(IntVector& fs) {
	Int count;
	IntVector Idf(mFacets.size(),0);
	IntVector Idc(mCells.size(),0);
	
	/*erase facet reference*/
	forEach(fs,i) {
		Int f = fs[i];
		Cell& co = mCells[mFOC[f]];
		Cell& coid = mFaceID[mFOC[f]];
		forEach(co,j) {
			if(co[j] == f) {
				co.erase(co.begin() + j); 
				coid.erase(coid.begin() + j);
				break; 
			}
		}
		Cell& cn = mCells[mFNC[f]];
		Cell& cnid = mFaceID[mFNC[f]];
		forEach(cn,j) {
			if(cn[j] == f) { 
				cn.erase(cn.begin() + j); 
				cnid.erase(cnid.begin() + j); 
				break; 
			}
		}
	}
	
	/*updated facet id*/
	forEach(fs,i)
		Idf[fs[i]] = Constants::MAX_INT;
	count = 0;
	forEach(mFacets,i) {
		if(Idf[i] != Constants::MAX_INT) 
			Idf[i] = count++;
		else
			mFacets[i].clear();
	}

	/*erase facets*/
	IntVector fzeroIndices;
	forEach(mFacets,i) {
		if(mFacets[i].size() == 0)
			fzeroIndices.push_back(i);
	}
	erase_indices(mFacets,fzeroIndices);
	erase_indices(mFOC,fzeroIndices);
	erase_indices(mFNC,fzeroIndices);
	erase_indices(mFC,fzeroIndices);
	erase_indices(mFN,fzeroIndices);
	/*updated cell id*/
	count = 0;
	forEach(mCells,i) {
		if(mCells[i].size() != 0) 
			Idc[i] = count++;
		else
			Idc[i] = Constants::MAX_INT;
	}
	/*erase cells*/
	IntVector czeroIndices;
	forEach(mCells,i) {
		if(mCells[i].size() == 0)
			czeroIndices.push_back(i);
	}
	erase_indices(mCells,czeroIndices);
	erase_indices(mFaceID,czeroIndices);
	erase_indices(mCC,czeroIndices);
	erase_indices(mCV,czeroIndices);
	
	/*updated facet id*/
	forEach(mCells,i) {
		forEach(mCells[i],j) {
			mCells[i][j] = Idf[mCells[i][j]];
		}
	}
	/*facet owner and neighbor*/
	forEach(mFacets,i) {
		mFOC[i] = Idc[mFOC[i]];
		mFNC[i] = Idc[mFNC[i]];
	}
	/*patches*/
	forEachIt(Boundaries,mBoundaries,it) {
		IntVector& mB = it->second;
		forEach(mB,i)
			mB[i] = Idf[mB[i]];
	}
}
/***************
* Refine Mesh
****************/

#define RDEBUG

void Mesh::MeshObject::straightEdges(const Facet& f, Facet& r, Facet& removed) {	
	Vector v1,v2,v3;
	Int fs = f.size();
	v1 = mVertices[f[fs - 1]];
	for(Int j = 0;j < fs;j++) {
		v2 = mVertices[f[j]];
		if(j + 1 == fs)
			v3 = mVertices[f[0]];
		else
			v3 = mVertices[f[j + 1]];
		Scalar m = mag((v2 - v1) ^ (v3 - v1));
		if(!equal(m,Scalar(0))) {
			r.push_back(f[j]);
			v1 = v2;
		} else
			removed.push_back(f[j]);
	}
}
bool Mesh::MeshObject::straightFaces(const Facet& f1_,const Facet& f2_) {
    Facet f1,f2,rr;
	straightEdges(f1_,f1,rr);
	straightEdges(f2_,f2,rr);
	Vector N1 = ((mVertices[f1[1]] - mVertices[f1[0]])
		^ (mVertices[f1[2]] - mVertices[f1[0]]));
	Vector N2 = ((mVertices[f2[1]] - mVertices[f2[0]])
		^ (mVertices[f2[2]] - mVertices[f2[0]]));
	Scalar e = mag(N1 ^ N2);
	if(equal(e,Scalar(0))) {
		Int count = 0;
		forEach(f1,i) {
			forEach(f2,j) {
				if(f1[i] == f2[j]) {
					count++;
					if(count == 2)
						return true;
				}
			}
		}
	}
	return false;
}
bool Mesh::MeshObject::mergeFacets(const Facet& f1_,const Facet& f2_, Facet& f) {
	//make them same direction
	Vector N1 = ((mVertices[f1_[1]] - mVertices[f1_[0]])
		^ (mVertices[f1_[2]] - mVertices[f1_[0]]));
	Vector N2 = ((mVertices[f2_[1]] - mVertices[f2_[0]])
		^ (mVertices[f2_[2]] - mVertices[f2_[0]]));
	Facet f1 = f1_;
	Facet f2 = f2_;
	if(dot(N1,N2) < 0) {
		forEach(f2_,j)
			f2[f2.size() - j - 1] = f2_[j];
	}
	//rotate vectors
	{
		Int v1 = -1,v2 = -1;
		forEach(f1,i) {
			Int j;
			for(j = 0;j < f2.size();j++) {
				if(f1[i] == f2[j]) 
					break;
			}
			if(j == f2.size()) {
				v1 = i;
				break;
			}
		}
		forEach(f2,i) {
			Int j;
			for(j = 0;j < f1.size();j++) {
				if(f2[i] == f1[j]) 
					break;
			}
			if(j == f1.size()) {
				v2 = i;
				break;
			}
		}
		std::rotate(f1.begin(),f1.begin() + v1,f1.end());
		std::rotate(f2.begin(),f2.begin() + v2,f2.end());
	}	
	//merge faces
	int a[2];
	int b[2];
	Int count;
	
	count = 0;
	forEach(f1,i) {
		forEach(f2,j) {
			if(f1[i] == f2[j]) {
				a[count] = i;
				b[count] = j;
				count++;
				if(count == 2) goto ENDM;
			}
		}
	}
	if(count < 2) 
		return false;
ENDM:;

    //merge
	f.clear();
	for(Int i = 0;i <= a[0];i++)
		f.push_back(f1[i]);
	for(Int i = b[0] + 1;i < f2.size();i++)
		f.push_back(f2[i]);
	for(Int i = 0;i < b[1];i++)
		f.push_back(f2[i]);
	for(Int i = a[1];i < f1.size();i++)
		f.push_back(f1[i]);

	//straighten edges
	Facet r,rr;
	straightEdges(f,r,rr);
	f = r;
	
	return true;
}
void Mesh::MeshObject::addVerticesToEdge(const int va, Facet& f, const Facet& fp) {

	Vector v1 = mVertices[f[f.size() - 1]];
	Vector v2 = mVertices[va];
	Scalar e1 = dot((v1 - v2),(v1 - v2));

	forEach(fp,i) {
		Vector v = mVertices[fp[i]];
		Scalar e = mag((v - v1)^(v2 - v1));
		if(equal(e,Scalar(0))) {
			e = dot((v - v2),(v1 - v2));
			if((e > Scalar(0)) && (e < e1)) {
				f.push_back(fp[i]);
			}
		}
	}
}
void Mesh::MeshObject::calcFaceCenter(const Facet& f,Vector& fCj) {
	{
		Vector C(0);
		forEach(f,j)
			C += mVertices[f[j]];
		fCj = C / Scalar(f.size());
	}
	Vector v1,v2,v3;
	Vector C(0);
	Scalar Ntot = Scalar(0);
	
	v1 = fCj;
	forEach(f,j) {
		v2 = mVertices[f[j]];
		if(j + 1 == f.size())
			v3 = mVertices[f[0]];
		else
			v3 = mVertices[f[j + 1]];
		Vector Ni = ((v2 - v1) ^ (v3 - v1));
		Scalar magN = mag(Ni);
		Vector Ci = magN * ((v1 + v2 + v3) / 3);


		C += Ci;
		Ntot += magN;
	}
	fCj = C / Ntot;    /*corrected face centre*/
}
void Mesh::MeshObject::calcCellCenter(const Cell& c, Vector& cCj) {
	
	//approximate cell centre
	Int cnt = 0;
	cCj = Vector(0);
	forEach(c,i) {
		Facet& f = mFacets[c[i]];
		forEach(f,j) {
			cCj += mVertices[f[j]];
			cnt++;
		}
	}
	cCj /= cnt;
	
	//exact cell centre
	Scalar Vt(0);
	Vector Ct(0);
	forEach(c,i) {
		Facet& f = mFacets[c[i]];
	
		Vector fCj,fNj,C;
	
		C = Vector(0);
		forEach(f,j)
			C += mVertices[f[j]];
		fCj = C / Scalar(f.size());
	
		Vector v1,v2,v3;
		Vector N(0);
		Scalar Ntot = Scalar(0);
	
		C = Vector(0);
		v1 = fCj;
		forEach(f,j) {
			v2 = mVertices[f[j]];
			if(j + 1 == f.size())
				v3 = mVertices[f[0]];
			else
				v3 = mVertices[f[j + 1]];
		
			Vector Ni = ((v2 - v1) ^ (v3 - v1));
			Scalar magN = mag(Ni);
			Vector Ci = magN * ((v1 + v2 + v3) / 3);
		
			C += Ci;
			Ntot += magN;
			N += Ni;
		}
		fCj = C / Ntot;
		fNj = N / Scalar(2);
	
		Scalar Vi = mag((cCj - fCj) & fNj);
		Ct += Vi * (3 * fCj + cCj) / 4;
		Vt += Vi;
	}

	cCj = Ct / Vt;
}
void Mesh::MeshObject::refineFacets(const IntVector& rFacets,IntVector& refineF, 
				 IntVector& startF, IntVector& endF,const Int ivBegin) {
					 
	using namespace Controls;
	
	const bool refine3D = equal(Scalar(0.0),mag(refine_params.dir));
	
	forEach(rFacets,i) {
		
		Int fi = rFacets[i];
		Facet f = mFacets[fi];
		Vector C;
		
		startF[fi] = mFacets.size();
		calcFaceCenter(f,C);
		
		/*add vertices*/
		mVertices.push_back(C);
		Int fci = mVertices.size() - 1;
		
		Facet fr;
		{
			Facet fs;
			straightEdges(f,fs,fr);
			f = fs;
		}

		IntVector midpts;
		midpts.resize(f.size(),Constants::MAX_INT);
		
		Vector v1,v2;
		Int fmid = Constants::MAX_INT;
		bool anisotropic = false;
		
		forEach(f,j) {
			v1 = mVertices[f[j]];
			Int nj = j + 1;
			if(nj == f.size())
				nj = 0;
			v2 = mVertices[f[nj]];

			if(!refine3D) {
				Scalar angle = acos(mag(dot(v2 - v1,refine_params.dir)) /
					          (mag(v2 - v1) * mag(refine_params.dir)));
				if(angle < Constants::PI/4 || angle >= 3 * Constants::PI/4) {
					anisotropic = true;
					continue;
				}
			}
			
			if(j < fmid) fmid = j;
			Vector Ce = (v1 + v2) / 2.0;

			//duplicate
			Int k = ivBegin;
			for(; k < mVertices.size();k++) {
				if(equal(Ce,mVertices[k])) {
					midpts[j] = k;
					break;
				}
			}
			if(k == mVertices.size() && fr.size()) {
				forEach(fr,i) {
					Int m = fr[i];
					if(equal(Ce,mVertices[m])) {
						midpts[j] = m;
						k = m;
						break;
					}
				}
			}
			if(k == mVertices.size()) {
				mVertices.push_back(Ce);
				midpts[j] = k;
			}
		}
		
		/*Add faces*/
#define ADD(x) {					\
	if(fr.size())					\
		addVerticesToEdge(x,fn,fr);	\
	fn.push_back(x);				\
}
		
		forEachS(f,j,fmid) {
			Facet fn;
			fn.push_back(midpts[j]);
			Int k = j + 1;
			for(;k < f.size();k++) {
				ADD(f[k]);
				Int k1 = midpts[k];
				if(k1 == Constants::MAX_INT);
				else { 
					ADD(k1);
					break;
				}
			}
			if(k == f.size()) {
				for(k = 0;k <= fmid;k++)
					ADD(f[k]);
				ADD(midpts[fmid]);
				j = f.size();
			} else {
				j = k - 1;
			}
			
			if(!anisotropic)
				fn.push_back(fci);

			mFacets.push_back(fn);
			startF.push_back(0);
			endF.push_back(0);
			refineF.push_back(0);
		}

#undef ADD
		
		endF[fi] = mFacets.size();
	}
}

void Mesh::MeshObject::refineCell(Cell& c,IntVector& cr,
				IntVector& refineF,IntVector& startF,IntVector& endF,
				Cells& newc,IntVector& delFacets
				) {

	Vector C;
	calcCellCenter(c,C);
	mVertices.push_back(C);
	Int cci = mVertices.size() - 1;
	Int ifBegin1 = mFacets.size();
	forEach(c,j) {
		Int fi = c[j];
	
		/*cell is refined but face is not?*/
		IntVector list;
		if(cr[j] != 1) {
			list.push_back(fi);
		} else {
			for(Int k = startF[fi]; k < endF[fi];k++)
				list.push_back(k);
		}

		/*refine cells*/
		forEach(list,k) {
			Int fni = list[k];
			Facet f = mFacets[fni];
	
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
				Int k = ifBegin1;
				for(; k < mFacets.size();k++) {
					if(equal(fn,mFacets[k])) {
						cn.push_back(k);
						break;
					}
				}
				if(k == mFacets.size()) {
					mFacets.push_back(fn);
					startF.push_back(0);
					endF.push_back(0);
					refineF.push_back(0);
					cn.push_back(k);
				}
			}
			newc.push_back(cn);
		}
	}
	/* *********************
	 *
	 * Merge cells 
	 *
	 * *********************/
	{
		/*cells on same face*/
		IntVector ownerf;
		forEach(c,j) {
			Int fi = c[j];
			if(cr[j] != 1) {
				ownerf.push_back(cr[j]);
			} else {
				for(Int k = startF[fi]; k < endF[fi];k++)
					ownerf.push_back(fi);
			}
#ifdef RDEBUG
			if(j == 0) std::cout << "====================================\n";
			std::cout << cr[j] << " ( " << startF[fi] << "  " << endF[fi] << " ) " << std::endl;
#endif
		}
	
		/* Find among new cells on different owner faces,
		 * and merge if they have a shared face
		 */
#ifdef RDEBUG
		std::cout << "============\n";
		std::cout << " Old cell " << c << "\ncC = " << C << std::endl;
		std::cout << "============\n";
		std::cout << " New cells before merge " << newc.size() << std::endl;
		std::cout << "============\n";
		std::cout << newc << std::endl;
#endif
		
		Cells mergec;
		forEach(newc,j) {
			Int o1 = ownerf[j];
			Cell& c1 = newc[j];
			IntVector mg;
			forEachS(newc,k,j+1) {
				Int o2 = ownerf[k];
				Cell& c2 = newc[k];
				if(o1 == o2) continue;
			
				/*Compare faces*/
				forEach(c1,m) {
					Int f1 = c1[m];
					forEach(c2,n) {
						Int f2 = c2[n];
						if(f1 == f2) {
							mg.push_back(k);
							goto END;
						}
					}
				}
				END:;
			}
			mergec.push_back(mg);
		}

		/*merge cells*/
		bool inserted;
		do {
			inserted = false;
			forEach(mergec,j) {
				IntVector& cm = mergec[j];
				forEach(cm,k) {
					forEachS(mergec,m,j+1) {
						IntVector& cn = mergec[m];
						bool has = false;
						forEach(cn,z) {
							if(cn[z] == cm[k]) {
								has = true;
								break;
							}
						}
						if(has) {
							inserted = true;
							cm.insert(cm.end(),cn.begin(),cn.end());
							cm.push_back(m);
							cn.clear();
						}
					}
				}
			}
		} while(inserted);
	
		IntVector erasei;
		erasei.assign(newc.size(),0);
		forEachR(mergec,j) {
			IntVector& cm = mergec[j];
			Cell& c1 = newc[j];
			forEach(cm,k) {
				if(erasei[cm[k]]) continue;
				erasei[cm[k]] = 1;
				Cell& c2 = newc[cm[k]];
			
				//remove internal faces
				forEach(c2,m) {
					Int c2m = c2[m];
					Int n;
					for(n = 0;n < c1.size();n++) {
						if(c2m == c1[n])
							break;
					}
					if(n == c1.size())
						c1.push_back(c2m);
					else {
						delFacets.push_back(c1[n]);
						c1.erase(c1.begin() + n);
					}
				}
			}
		}
		/*erase cells*/
		IntVector erasec;
		forEach(erasei,j) {
			if(erasei[j])
				erasec.push_back(j);
		}
		erase_indices(newc,erasec);
	
#ifdef RDEBUG
		std::cout << " New cells " << newc.size() << std::endl;
		std::cout << "============\n";
		std::cout << newc << std::endl;
#endif
	
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
				//more than one face shared between two cells
				if(shared1.size() > 1) {
				
					Facet f = mFacets[c1[shared1[0]]];
					IntVector flag;
					flag.assign(shared1.size(),0);
					bool has,mhas,bmerge = false;
					do {
						has = false;
						mhas = false;
						forEachS(shared1,p,1) {
							if(flag[p]) continue;
							Facet& f1 = mFacets[c1[shared1[p]]];
							Facet f2 = f;
							if(mergeFacets(f2,f1,f)) {
								flag[p] = 1;
								mhas = true;
								bmerge = true;
							} else {
								has = true;
							}
						}
					} while(has && mhas);
#ifdef RDEBUG
					if(f.size() > 4)  {
						std::cout << "Merge failed.\n";
						std::cout << "=======================\n";
						{
							Facet r,rr;
							straightEdges(f,r,rr);
							std::cout << "Merged face  : " <<  f << std::endl;
							std::cout << "Straight Edge: " << r << std::endl;
							std::cout << "Removed      : " << rr << std::endl;
						}
						forEach(f,i)
							std::cout << mVertices[f[i]] << " : ";
						std::cout << std::endl;
						
						std::cout << newc.size() << " " << shared1 << " @@@@@ " << shared2 << std::endl;
						std::cout << c1 << std::endl;
						std::cout << c2 << std::endl;
						forEach(shared1,p) {
							Facet& f1 = mFacets[c1[shared1[p]]];
							forEach(f1,i)
								std::cout << " ( " << mVertices[f1[i]] << " ) ";
							std::cout << " = " << f1 << std::endl;
						}
						std::cout << "=======================\n";
						exit(0);
					}
#endif
					//remove merged facets
					forEach(shared1,p)
						delFacets.push_back(c1[shared1[p]]);
					sort(shared2.begin(),shared2.end());
					erase_indices(c1,shared1);
					erase_indices(c2,shared2);
					
					//add new face
					Int m = ifBegin1;
					for(; m < mFacets.size();m++) {
						if(equal(f,mFacets[m]))
							break;
					}
					if(m == mFacets.size()) {
						mFacets.push_back(f);
						startF.push_back(0);
						endF.push_back(0);
						refineF.push_back(0);
						m = mFacets.size() - 1;
					} else {
						std::cout << "****** Found ******\n";
					}
					c1.push_back(m);
					c2.push_back(m);
				}
			}
		}
	}
#ifdef RDEBUG
	std::cout << " After merge new cells " << newc.size() << std::endl;
	std::cout << "============================\n";
	std::cout << newc << std::endl;
	forEach(newc,i) {
		Cell& c = newc[i];
		forEach(c,j) {
			if(mFacets[c[j]].size() <= 3)
				exit(0);
		}
	}
#endif

}

void Mesh::MeshObject::refineMesh(IntVector& rCells,IntVector& cellMap) {
	using namespace Controls;
	
#ifdef RDEBUG
	std::cout << rCells << std::endl;
	std::cout << "----------\n";
#endif
	
	/*Initialize coarse to fine grid map*/
	cellMap.assign(mBCS,0);
	forEach(cellMap,i)
		cellMap[i] = i;
	
	/*Initialize face refinement info*/
	IntVector refineF;
	Cells crefineF;
	
	refineF.assign(mFacets.size(),0);
	
	forEach(rCells,i) {
		IntVector v = mCells[rCells[i]];
		forEach(v,k) v[k] = 0;
		crefineF.push_back(v);
	}

	forEach(rCells,i) {
		Int ci = rCells[i];
		Cell& c = mCells[ci];
		Cell& cr = crefineF[i];
		forEach(c,j) {
			Int fj = c[j];
			bool straight = false;
			forEach(c,k) {
				if(j == k) continue;
				bool b1 = straightFaces(mFacets[c[j]],mFacets[c[k]]);
				if(b1) {
					Int nid = ((j > k) ? j : k) + 2;
					if(nid > cr[j]) 
						cr[j] = nid;
					straight = true;
				}
			}
			if(!straight) {
				refineF[fj] = 1;
				cr[j] = 1;
			}
		}
		forEach(cr,j) {
			if(cr[j] >= 2)
				cr[j] = cr[cr[j] - 2];
		}
	
		const Int HALF = Constants::MAX_INT / 2;
		IntVector flag;
		flag.assign(cr.size(),0);
		forEach(cr,k) {
			if(flag[k]) continue;
		
			Int crk = cr[k];
			if(crk < 2) continue;
			Vector C(0);
			Int cnt = 0;
			forEach(cr,j) {
				Int crj = cr[j];
				if(crj >= HALF)
					crj -= HALF;
				if(crk == crj) {
					C += mFC[c[j]];
					cnt++;
				}
			}
			C /= cnt;
			if(cnt <= 2) continue;

			Scalar e = dot(mFC[c[k]] - C, mFC[c[crk - 2]] - C);
			if(e > Scalar(0))
				cr[k] += HALF;
		
			Scalar maxd = -1e30,mind = 1e30;
			Int ix,in;
			forEach(cr,j) {
				Int crj = cr[j];
				if(crk == crj) {
					flag[j] = 1;
					Scalar e = dot(mFC[c[j]] - C, mFC[c[crk - 2]] - C);
					if(e > Scalar(0)) {
						if(e < mind) mind = e;
						in = j;
					} else {
						if(e > maxd) maxd = e;
						ix = j;
					}
				}
			}

			forEachS(cr,j,k) {
				Int crj = cr[j];
				if(crk == crj)
					cr[j] = crj + c[j];
			}
		
			cr[in] = cr[ix];
		}
	}
		
	/***************************
	 * Refine facets
	 ***************************/
	Int ivBegin = mVertices.size();
	IntVector startF;
	IntVector endF;
	startF.assign(mFacets.size(),0);
	endF.assign(mFacets.size(),0);
	
	IntVector rFacets;
	{
		forEach(refineF,i) {
			if(refineF[i])
				rFacets.push_back(i);
		}
		
		refineFacets(rFacets,refineF,startF,endF,ivBegin);
	}
	/*************************
	 * Refine cells
	 *************************/
	IntVector delFacets;
	
	forEach(rCells,i) {
		Cells newc;
		Cell cr;
		Vector C;

		Int ci = rCells[i];
		newc.push_back(mCells[ci]);
		cr = crefineF[i];
		C = gcC[ci];

		const int NR = 1;
		for(Int m = 0;m < NR;m++) {
			
#ifdef RDEBUG
			if(m == 0) {
				std::cout << "========================================================================\n";
				std::cout << "Cell " << i << std::endl;
				std::cout << "==================\n";
				std::cout << newc << std::endl;
			}
			std::cout << "Refinement level " << m + 1 << std::endl;
			std::cout << "==================\n";
#endif
			
			if(m) {
				IntVector rfFacets;
				forEach(newc,j)
					rfFacets.insert(rfFacets.end(),newc[j].begin(),newc[j].end());
				std::sort(rfFacets.begin(), rfFacets.end());
				rfFacets.erase(std::unique(rfFacets.begin(),rfFacets.end()), rfFacets.end());
				
				IntVector alRef;
				forEach(rfFacets,j) {
					Int fi = rfFacets[j];
					if(refineF[fi])
						alRef.push_back(j);
					else {
						refineF[fi] = 1;
						rFacets.push_back(fi);
					}
				}
				erase_indices(rfFacets,alRef);

				refineFacets(rfFacets,refineF,startF,endF,ivBegin);
			}

			Cells newcm;
			forEach(newc,b) {
				Cell& c = newc[b];
				if(m) cr.assign(c.size(),1);
				
				Cells newcn;
				refineCell(c,cr,refineF,startF,endF,newcn,delFacets);
				newcm.insert(newcm.end(),newcn.begin(),newcn.end());
			}

			newc = newcm;
			
#ifdef RDEBUG
			std::cout << newc << std::endl;
#endif
		}
		
		/*add cells*/
		mCells.insert(mCells.end(),newc.begin(),newc.end());
		forEach(newc,j)
			cellMap.push_back(ci);
	}

	/*************************************
	 *  Remove refined facets and cells
	 *************************************/

	/*add additional facets to remove*/
	rFacets.insert(rFacets.end(),delFacets.begin(),delFacets.end());
	forEach(delFacets,i) {
		Int fi = delFacets[i];
		refineF[fi] = Constants::MAX_INT - 1;
	}
	
	/*adjust facet indexes in cells and boundaries*/
	Int count = 0;
	forEach(refineF,i) {
		if(refineF[i] == 0) refineF[i] = count++;
		else if(refineF[i] == Constants::MAX_INT - 1);
		else refineF[i] = Constants::MAX_INT;
	}

	/*erase old facets and add the new ones*/
#define ERASEADD() {										\
	bool has;												\
	do {													\
		has = false;										\
		IntVector newF,eraseF;								\
		forEach(c,j) {										\
			Int fi = c[j];									\
			if(refineF[fi] == Constants::MAX_INT) {			\
				eraseF.push_back(j);						\
				for(Int k = startF[fi]; k < endF[fi];k++) {	\
					newF.push_back(k);						\
					has = true;								\
				}											\
			}												\
		}													\
		erase_indices(c,eraseF);							\
		c.insert(c.end(),newF.begin(),newF.end());			\
	} while (has);											\
															\
	forEach(c,j) 											\
		c[j] = refineF[c[j]];								\
}
	forEach(mCells,i) {
		Cell& c = mCells[i];
		ERASEADD();
	}
	forEachIt(Boundaries,mBoundaries,it) {
		IntVector& c = it->second;
		ERASEADD();
	}
#undef ERASEADD
	
	/*erase facets*/
	sort(rFacets.begin(),rFacets.end());
	erase_indices(mFacets,rFacets);
	
	/*erase cells*/
	erase_indices(mCells,rCells);
	erase_indices(cellMap,rCells);
	mBCS = mCells.size();
	
	/*******************************************************
	 * a) Remove unused vertices
	 * b) Break edges of faces that are not set for refinement 
	 *    but should be due to neighboring refined faces 
	 *******************************************************/
	{
		/* remove unused vertices */
		IntVector isUsed(mVertices.size(),0);
		forEach(mFacets,i) {
			Facet& f = mFacets[i];
			forEach(f,j)
				isUsed[f[j]] = 1;
		}
		IntVector rVertices;
		Int cnt = 0;
		forEach(isUsed,i) {
			if(!isUsed[i])
				rVertices.push_back(i);
			else
				isUsed[i] = cnt++;
		}
		if(rVertices.size()) {
			forEach(mFacets,i) {
				Facet& f = mFacets[i];
				forEach(f,j)
					f[j] = isUsed[f[j]];
			}
			erase_indices(mVertices,rVertices);
		}
		
		/* break edges of faces */
		bool hasMid = false;
		forEach(mFacets,i) {
			Facet& f = mFacets[i];
			IntVector addf;
			Vector v1,v2;
			addf.assign(f.size(),Constants::MAX_INT);
			forEach(f,j) {
				v1 = mVertices[f[j]];
				Int nj = j + 1;
				if(nj == f.size())
					nj = 0;
				v2 = mVertices[f[nj]];
				Vector Ce = (v1 + v2) / 2.0;

				//duplicate
				Int k = ivBegin;
				for(; k < mVertices.size();k++) {
					if(equal(Ce,mVertices[k])) {
						hasMid = true;
						addf[j] = k;
						break;
					}
				}
			}
			if(hasMid) {
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
	/****************
	 *   End
	 ****************/
}

