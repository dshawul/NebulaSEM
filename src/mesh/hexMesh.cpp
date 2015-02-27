#include "hexMesh.h"
#include <cmath>

using namespace Mesh;

/**
Center of a triangle
*/
Vector center(const Vector& v1,const Vector& v2,const Vector& v3) {
	Vector v12 = v1 - v2;
	Vector v13 = v1 - v3;
	Vector v23 = v2 - v3;
	Scalar d = 2 * magSq(v12 ^ v23);
	Scalar a = magSq(v23) * (v12 & v13) / d;
	Scalar b = magSq(v13) * (-v12 & v23) / d;
	Scalar c = magSq(v12) * (v13 & v23) / d;
	return a * v1 + b * v2 + c * v3;
}
/**
Add different shapes of edges
*/
void ADDV(int w,Scalar m,Edge* edges,Vector* vd) {
	Edge& e = edges[w];
	if(e.type == NONE) {
		vd[w] = (1 - m) * e.v[0] + (m) * e.v[1];
	} else if(e.type == ARC) {
		vd[w] = rotate(e.v[0] - e.v[3],e.N,e.theta * m) + e.v[3];
	} else if(e.type == COSINE) {
		vd[w] = (1 - m) * e.v[0] + (m) * e.v[1] + 
			pow(cos(3.1416 * (m - 0.5)),2) * e.N;
	} else if(e.type == QUAD) {
		vd[w] = (1 - m) * e.v[0] + (m) * e.v[1] + 
			(4 * m * (1 - m)) * e.N;
	}
}
/**
Generate hexahedral mesh
*/
void hexMesh(Int* n,Scalar* s,Int* type,Vector* vp,Edge* edges,MeshObject& mo) {
	Int i,j,k,m;

	/*for wall division set twice 
	number of divisions requested*/
	for(j = 0;j < 3;j++) {
		bool found = false;
		for(i = j;i < 12;i+=3) {
			if(type[i] == WALL) {
				if((n[j] % 2) && (n[j] != 1)) {
					found = true;
					break;
				}
			}
		}
		if(found) {
			n[j]++;
			for(i = j;i < 12;i+=3) {
				s[i] = 1 / s[i];
			}
		}
	}
	
	/*calculate scale*/
	Scalar* sc[12];
	for(i = 0;i < 12;i++) {
		Int nt = n[i / 4];
		sc[i] = new Scalar[nt + 1];
		if(type[i] == WALL)
			s[i] = pow(s[i],Scalar(1./(nt / 2.)));
		else
			s[i] = pow(s[i],Scalar(1./nt));
	}
	for(i = 0;i < 12;i++) {
		Int nt = n[i / 4];
		Scalar r = s[i];
		if(nt == 1) {
			sc[i][0] = 0;
			sc[i][1] = 1;
		} else {
			if(type[i] == WALL)
				nt /= 2;
			for(j = 0;j <= nt;j++) {
				if(equal(r,Scalar(1))) 
					sc[i][j] = Scalar(j) / (nt);
				else
					sc[i][j] = (1 - pow(r,Scalar(j))) / (1 - pow(s[i],Scalar(nt)));
			}
			if(type[i] == WALL) {
				for(j = 0;j <= nt;j++)
					sc[i][j] /= 2;
				for(j = 0;j <= nt;j++)
					sc[i][j + nt] = Scalar(1.0) - sc[i][nt - j];
			}
		}
	}
	for(i = 0;i < 12;i++) {
		Edge& e = edges[i];
		if(e.type == ARC) {
			Vector C = center(e.v[0],e.v[1],e.v[2]);
			Vector r1 = e.v[0] - C;
			Vector r2 = e.v[1] - C;
			e.theta = acos((r1 & r2) / (mag(r1) * mag(r2)));
			e.v[3] = C;
			e.N = (e.v[2] - e.v[0]) ^ (e.v[1] - e.v[0]);
			e.N = unit(e.N);
		} else if(e.type == COSINE || e.type == QUAD) {
			Vector mid = (e.v[1] + e.v[0]) / 2;
			e.N = e.v[2] - mid;
			e.L = mag(mid - e.v[0]) / 2;
		}
	}
    /*variables*/
	Int nx = n[0] + 1 , ny = n[1] + 1 , nz = n[2] + 1;
	const Int B1 = (nx - 0) * (ny - 1) * (nz - 1);
	const Int B2 = (nx - 1) * (ny - 0) * (nz - 1);
	const Int B3 = (nx - 1) * (ny - 1) * (nz - 0);
	IntVector VI(nx * ny * nz,0);
	IntVector FI(B1 + B2 + B3,0);

	/*vertices*/
	Vertex v,v1,v2,vd[12],vf[6];
	Scalar rx,ry,rz;

#define I0(i,j,k)  (i * ny * nz + j * nz + k)

#define ADDF(w,rr,rs,i00,i01,i10,i11,ir0,ir1,i0s,i1s) {	\
	vf[w] = Interpolate_face(							\
			rr,rs,										\
			vp[i00],vp[i01],vp[i10],vp[i11],			\
			vd[ir0],vd[ir1],vd[i0s],vd[i1s]);			\
}

#define ADDC() {										\
	v = Interpolate_cell(								\
			rx,ry,rz,									\
			vp[0],vp[4],vp[3],vp[7],					\
			vp[1],vp[5],vp[2],vp[6],					\
			vd[0],vd[3],vd[1],vd[2],					\
			vd[4],vd[7],vd[5],vd[6],					\
			vd[8],vd[11],vd[9],vd[10],					\
			vf[4],vf[5],vf[2],vf[3],vf[0],vf[1]);		\
}

#define ADD() {                                     \
	ADDV(0,sc[0][i],edges,vd);						\
	ADDV(1,sc[1][i],edges,vd);						\
	ADDV(2,sc[2][i],edges,vd);						\
	ADDV(3,sc[3][i],edges,vd);						\
	ADDV(4,sc[4][j],edges,vd);						\
	ADDV(5,sc[5][j],edges,vd);						\
	ADDV(6,sc[6][j],edges,vd);						\
	ADDV(7,sc[7][j],edges,vd);						\
	ADDV(8,sc[8][k],edges,vd);						\
	ADDV(9,sc[9][k],edges,vd);						\
	ADDV(10,sc[10][k],edges,vd);					\
	ADDV(11,sc[11][k],edges,vd);					\
	rx = i / Scalar(nx - 1);						\
	ry = j / Scalar(ny - 1);						\
	rz = k / Scalar(nz - 1);						\
	ADDF(0, rx,ry, 0,3,1,2, 0,1,4,5);				\
	ADDF(1, rx,ry, 4,7,5,6, 3,2,7,6);				\
	ADDF(2, rx,rz, 0,4,1,5, 0,3,8,9);				\
	ADDF(3, rx,rz, 3,7,2,6, 1,2,11,10);				\
	ADDF(4, ry,rz, 0,4,3,7, 4,7,8,11);				\
	ADDF(5, ry,rz, 1,5,2,6, 5,6,9,10);				\
	ADDC();											\
};

	/*interior*/
	for(j = 1;j < ny - 1;j++) {
		for(i = 1;i < nx - 1;i++) {
			for(k = 1;k < nz - 1;k++) {
				ADD();
				mo.v.push_back(v);
				VI[I0(i,j,k)] = mo.v.size() - 1;
			}
		}
	}
	mo.nv = mo.v.size();

	/*boundaries*/
	for(i = 0;i < nx; i += (nx - 1)) {
		for(j = 0;j < ny;j++) {
			for(k = 0;k < nz;k++) {
				ADD();
				mo.v.push_back(v);
				VI[I0(i,j,k)] = mo.v.size() - 1;
			}
		}
	}
	for(j = 0;j < ny; j += (ny - 1)) {
		for(i = 1;i < nx - 1;i++) {
			for(k = 0;k < nz;k++) {
				ADD();
				mo.v.push_back(v);
				VI[I0(i,j,k)] = mo.v.size() - 1;
			}
		}
	}
	for(k = 0;k < nz; k += (nz - 1)) {
		for(i = 1;i < nx - 1;i++) {
			for(j = 1;j < ny - 1;j++) {
				ADD();
				mo.v.push_back(v);
				VI[I0(i,j,k)] = mo.v.size() - 1;
			}
		}
	}
	/*end*/
#undef ADD
#undef ADDF
#undef ADDE

	delete[] sc[0];
	delete[] sc[1];
	delete[] sc[2];

	/*faces*/
#define I1(i,j,k) 	(i * (ny - 1) * (nz - 1) + j * (nz - 1) + k)
#define I2(i,j,k) 	(i * (ny - 0) * (nz - 1) + j * (nz - 1) + k + B1)
#define I3(i,j,k) 	(i * (ny - 1) * (nz - 0) + j * (nz - 0) + k + B1 + B2)

#define ADD(a1,a2,a3,a4) {							\
	Facet f;										\
	m = I0(i,j,k);									\
	f.push_back(VI[a1]);							\
	f.push_back(VI[a2]);							\
	f.push_back(VI[a3]);							\
	f.push_back(VI[a4]);							\
	mo.f.push_back(f);								\
};

	/*interior*/
	for(i = 0;i < nx - 1;i++) {
		for(j = 0;j < ny - 1;j++) {
			for(k = 1;k < nz - 1;k++) {
				ADD(m, m + ny * nz,m + ny * nz + nz, m + nz);
				FI[I3(i,j,k)] = mo.f.size() - 1;
			}
		}
	}
	for(i = 0;i < nx - 1;i++) {
		for(j = 1;j < ny - 1;j++) {
			for(k = 0;k < nz - 1;k++) {
				ADD(m,m + ny * nz,m + ny * nz + 1,m + 1);
				FI[I2(i,j,k)] = mo.f.size() - 1;
			}
		}
	}
	for(i = 1;i < nx - 1;i++) {
		for(j = 0;j < ny - 1;j++) {
			for(k = 0;k < nz - 1;k++) {
				ADD(m,m + nz,m + nz + 1,m + 1);
				FI[I1(i,j,k)] = mo.f.size() - 1;
			}
		}
	}
	mo.nf = mo.f.size();
	/*boundaries*/
	for(k = 0;k < nz; k += (nz - 1)) {
		for(i = 0;i < nx - 1;i++) {
			for(j = 0;j < ny - 1;j++) {
				ADD(m, m + ny * nz,m + ny * nz + nz, m + nz);
				FI[I3(i,j,k)] = mo.f.size() - 1;
			}
		}
	}
	for(j = 0;j < ny;j += (ny - 1)) {
		for(i = 0;i < nx - 1;i++) {
			for(k = 0;k < nz - 1;k++) {
				ADD(m,m + ny * nz,m + ny * nz + 1,m + 1);
				FI[I2(i,j,k)] = mo.f.size() - 1;
			}
		}
	}
	for(i = 0;i < nx; i += (nx - 1)) {
		for(j = 0;j < ny - 1;j++) {
			for(k = 0;k < nz - 1;k++) {
				ADD(m,m + nz,m + nz + 1,m + 1);
				FI[I1(i,j,k)] = mo.f.size() - 1;
			}
		}
	}
	
	/*end*/
#undef ADD

	/*cells*/
	for(i = 0;i < nx - 1;i++) {
		for(j = 0;j < ny - 1;j++) {
			for(k = 0;k < nz - 1;k++) {
				Cell c;
				m = I3(i,j,k);
				c.push_back(FI[m]);
				c.push_back(FI[m + 1]);
				
				m = I2(i,j,k);
				c.push_back(FI[m]);
				c.push_back(FI[m + (nz - 1)]);
				
				m = I1(i,j,k);
				c.push_back(FI[m]);
				c.push_back(FI[m + (ny - 1) * (nz - 1)]);

				mo.c.push_back(c);
			}
		}
	}
	mo.nc = mo.c.size();
#undef I0
#undef I1
#undef I2
#undef I3
	/*remove duplicates*/
	int deformed = 0;
	for(i = 0;i < 8;i++) {
		for(j = i + 1;j < 8;j++) {
			if(equal(vp[i],vp[j])) {
				deformed = 1;
				break;
			}
		}
	}
	if(deformed)
		remove_duplicate(mo);
	/*end*/
}

/**
Remove duplicate vertices,faces and cells
*/
void remove_duplicate(Mesh::MeshObject& p) {
	Int i,j,sz,corr;
	int count;
	/*vertices*/
	sz = p.v.size();
	corr = 0;
	std::vector<int> dup(sz,0);
	for(i = 0;i < sz;i++) {
		for(j = sz - 1;j >= i + 1;j--) {
			if(equal(p.v[i],p.v[j])) {
				dup[i] = -int(j);
				if(i < p.nv) corr++;
				break;
			}
		}
	}
	p.nv -= corr;
	//remove duplicate vertices
	{
		Vertices vt(p.v.begin(), p.v.end());
		p.v.clear();
		count = 0;
		for(i = 0;i < sz;i++) {
			if(!dup[i]) {
				p.v.push_back(vt[i]);
			    dup[i] = count++;
			}
		}
		for(i = 0;i < sz;i++) {
			if(dup[i] < 0) 
				dup[i] = dup[-dup[i]];
		}
	}
	/*faces*/
	sz = p.f.size();
	for(i = 0;i < sz;i++) {
		Facet& f = p.f[i];
		forEach(f,j)
			f[j] = dup[f[j]];
	}
	dup.clear();
	dup.assign(sz,0);
	count = 0;
	corr = 0;
	for(i = 0;i < sz;i++) {
		Facet& f = p.f[i];
		forEach(f,j) {
			forEachS(f,k,j+1) {
				if(f[j] == f[k]) {
					f.erase(f.begin() + k);
					k--;
				}
			}
		}
		if(f.size() < 3) {
			dup[i] = -1;
			if(i < p.nf) corr++;
		} else {
			dup[i] = count;
			count++;
		}
	}
	p.nf -= corr;
	//remove deformed faces
	{
		Facets ft(p.f.begin(), p.f.end());
		p.f.clear();
		for(i = 0;i < sz;i++) {
			if(dup[i] >= 0) p.f.push_back(ft[i]);
		}
	}
	/*cells*/
	sz = p.c.size();
	for(i = 0;i < sz;i++) {
		Cell& c = p.c[i];
		forEach(c,j) {
			if(dup[c[j]] < 0) {
				c.erase(c.begin() + j);
				j--;
			} else
				c[j] = dup[c[j]];
		}
	}
}

#define MAXNUM 1073741824

/**
Merge mesh m2 onto m1 (internal) and b (boundary) meshes
*/
void merge(MeshObject& m1,MergeObject& b,MeshObject& m2) {
	Int i,j,found,s0,s1,s2,s3;

    //vertices
	{
		s0 = m1.v.size();
		s1 = m2.nv;
		s2 = m2.v.size();
		s3 = b.vb.size();
		m1.v.insert(m1.v.end(),m2.v.begin(),m2.v.begin() + s1);

		IntVector locv(s2 - s1,MAXNUM);
		for(i = s1;i < s2;i++) {
			found = 0;
			for(j = 0;j < s3;j++) {
				if(equal(m2.v[i],b.vb[j])) {
					locv[i - s1] += j;
					found = 1;
					break;
				}
			}
			if(!found) {
				b.vb.push_back(m2.v[i]);
				locv[i - s1] += b.vb.size() - 1;
			}
		}
		forEach(m2.f,i) {
			Facet& ft = m2.f[i];
			forEach(ft,j) {
				if(ft[j] >= s1) {
					ft[j] = locv[ft[j] - s1];
				} else {
					ft[j] += s0;
				}
			}
		}
	}
	//faces
	{
		s0 = m1.f.size();
		s1 = m2.nf;
		s2 = m2.f.size();
		s3 = b.fb.size();
		m1.f.insert(m1.f.end(),m2.f.begin(),m2.f.begin() + s1);

		IntVector index0(s3,0),index1(s2 - s1,0);
		Int count = 0;
		b.fb.reserve(s3 + s2 - s1);
		for(j = 0;j < s3;j++) {
			found = 0;
			for(i = s1;i < s2;i++) {
				if(!index1[i - s1] && equal(m2.f[i],b.fb[j])) {

					m1.f.push_back(b.fb[j]);
					index0[j]      = m1.f.size() - 1;
					index1[i - s1] = m1.f.size() - 1;

					found = 1;
					break;
				}
			}
			if(!found) {
				index0[j] = MAXNUM + count;
				b.fb[count] = b.fb[j];
				count++;
			}
		}
		for(i = s1;i < s2;i++) {
			if(!index1[i - s1]) {
				index1[i - s1] = MAXNUM + count;

				if(count >= s3) b.fb.push_back(m2.f[i]);
				else b.fb[count] = m2.f[i];
				count++;
			}
		}
		b.fb.resize(count);

		forEach(m1.c,i) {
			Cell& ct = m1.c[i];
			forEach(ct,j) {
				if(ct[j] >= MAXNUM) {
					ct[j] = index0[ct[j] - MAXNUM];
				}
			}
		}
		forEach(m2.c,i) {
			Cell& ct = m2.c[i];
			forEach(ct,j) {
				if(ct[j] >= s1) {
					ct[j] = index1[ct[j] - s1];
				} else {
					ct[j] += s0;
				}
			}
		}
	}
	//cells
	{
		m1.c.insert(m1.c.end(),m2.c.begin(),m2.c.end());
	}
}
/**
Merge boundary and internals
*/
void merge(Mesh::MeshObject& m,MergeObject& b) {
	m.nv = m.v.size();
	m.nf = m.f.size();
	m.nc = m.c.size();

	m.v.insert(m.v.end(),b.vb.begin(),b.vb.end());
	m.f.insert(m.f.end(),b.fb.begin(),b.fb.end());
	forEach(m.f,i) {
		Facet& ft = m.f[i];
		forEach(ft,j) {
			if(ft[j] >= MAXNUM) {
				ft[j] -= MAXNUM;
				ft[j] += m.nv;
			}
		}
	}
	forEach(m.c,i) {
		Cell& ct = m.c[i];
		forEach(ct,j) {
			if(ct[j] >= MAXNUM) {
				ct[j] -= MAXNUM;
				ct[j] += m.nf;
			}
		}
	}
}

#undef MAXNUM
