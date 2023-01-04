#include "field.h"
#include "system.h"
#include "metis.h"

using namespace std;

/**
  Define some geometric fields
 */
namespace Mesh {
    VectorVertexField vC(false);
    VectorFacetField  fC(false);
    VectorCellField   cC(false);
    VectorFacetField  fN(false);
    ScalarCellField   cV(false);
    ScalarFacetField  fI(false);
    ScalarFacetField  fD(false);
    ScalarCellField   yWall(false);
    IntVector         FO;
    IntVector         FN;
    IntVector  probeCells;
    Int         gBCSfield;
    Int         gBCSIfield;
    Int         gBFSfield;

    Vertices                 probePoints;
    vector<BasicBCondition*> AllBConditions;
}

std::list<BaseField*> BaseField::allFields;
std::vector<std::string> BaseField::fieldNames;

/**
  Refinement and domain decomposition parameters
 */
namespace Controls {
    RefineParams refine_params;
    DecomposeParams decompose_params;
}

/**
  Solver control parameters
 */
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
    FILE_FORMAT write_format = BINARY;
}
/**
  Find the last refined grid
 */
static int findLastRefinedGrid(Int step_) {
    int step = step_;
    for(;step >= 0;--step) {
        stringstream path;
        path << Mesh::gMeshName << "_" << step;
        string str = path.str();
        if(System::exists(str+".txt") ||
           System::exists(str+".bin"))
            break;
    }
    if(step < 0) step = step_;
    return step;
}
/**
  Load mesh
 */
bool Mesh::LoadMesh(Int step_, bool remove_empty) {
    /*load refined mesh*/
    int step = findLastRefinedGrid(step_);

    /*open grid file*/
    stringstream path;
    path << gMesh.name << "_" << step;
    string str = path.str();

    bool loaded = false, found = false;
    Int iter = 0;
    do {
        if(System::exists(str + ".txt")) {
            std::ifstream is(str + ".txt");
            loaded = gMesh.readTextMesh(is);
            found = true;
        } else if(System::exists(str + ".bin")) {
            Util::ifstream_bin is(str + ".bin");
            loaded = gMesh.readTextMesh(is);
            found = true;
        }
        path.clear();
        path << gMesh.name << "_" << 0;
        str = path.str();
        iter++;
    } while(found == false && iter == 1);

    /*load mesh*/
    if(loaded) {
        /*clear bc and probing points list*/
        Mesh::clearBC();
        Mesh::probePoints.clear();
        /*print info*/
        if(MP::printOn)
            cout << "--------------------------------------------\n";
        MP::printH("\t%d vertices\t%d facets\t%d cells\n",
                gVertices.size(),gFacets.size(),gCells.size());
        /*initialize mesh*/
        gMesh.addBoundaryCells();
        gMesh.calcGeometry();
        DG::init_poly();
        /* remove empty faces*/
        if(remove_empty) {
            auto it = gBoundaries.find("delete");
            if(it != gBoundaries.end()) {
                IntVector& fs = gBoundaries["delete"];
                gMesh.removeBoundary(fs);
                gBoundaries.erase(it);
                gMesh.removeUnusedVertices();
            }
        }
        /*erase interior and empty boundaries*/
        for(auto it = gBoundaries.begin();it != gBoundaries.end();) {
            if(it->second.size() <= 0 || 
                    it->first.find("interior") != std::string::npos
              ) {
                gBoundaries.erase(it++);
            } else ++it;
        }
        /*geometric mesh fields*/
        remove_fields();
        initGeomMeshFields();
        /*Adjust interpolation factor for "delete" faces*/
        if(!remove_empty) {
            auto it = gBoundaries.find("delete");
            if(it != gBoundaries.end()) {
                IntVector& fs = gBoundaries["delete"];
                forEach(fs,i) {
                    fI[fs[i]] = 1;
                }
            }
        }
        if(MP::printOn) 
            cout << "--------------------------------------------\n";
        return true;
    }
    return false;
}
/**
  Initialize geometric mesh fields
 */
void Mesh::initGeomMeshFields() {
    FO = gFOC;
    FN = gFNC;
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
    forEach(gCC,i) {
        cC[i] = gCC[i];
        cV[i] = gCV[i];
    }
    forEach(gFC,i) {
        fC[i] = gFC[i];
        fN[i] = gFN[i];
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
    /*finish comm*/
    commv.recv();
    commc.recv();
    /* Facet interpolation factor to the owner of the face.
     * Neighbor takes (1 - f) */
    forEach(gFacets,faceid) {
        Vector v0 = gVertices[gFacets[faceid][0]];
        for(Int n = 0; n < DG::NPF;n++) {
            Int k = faceid * DG::NPF + n;
            Int c1 = FO[k];
            Int c2 = FN[k];
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
    /*construct diffusivity factor*/
    if(DG::NPMAT) {
        using namespace DG;

        //penalty
        Scalar k = (NPX > NPY) ? NPX : 
            ((NPY > NPZ) ? NPY : NPZ);
        Scalar num;

        if(NP == k)
            num = (k + 1) * (k + 1) / 1;
        else if(NPX == 1 || NPY == 1 || NPZ == 1)
            num = (k + 1) * (k + 2) / 2;
        else
            num = (k + 1) * (k + 3) / 3;

        fD = cds(cV);
        forEach(fN,i)
            fD[i] = -num * mag(fN[i]) / (fD[i]);

        //diffusivity
        VectorCellField grad_psi = Vector(0);

        for(Int ci = 0; ci < gBCS; ci++) {
            forEachLglBound(ii,jj,kk) {
                Int index = INDEX4(ci,ii,jj,kk);

#define PSID(im,jm,km) {                    \
    Int index1 = INDEX4(ci,im,jm,km);       \
    Vector dpsi_ij;                         \
    DPSIR(dpsi_ij,im,jm,km);                \
    dpsi_ij = dot(Jinv[index1],dpsi_ij);    \
    grad_psi[index] += dpsi_ij;             \
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
            Int c1 = FO[i];
            Int c2 = FN[i];
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
        forEachIt(gBoundaries,it) {
            string bname = it->first;
            bc = new BCondition<Scalar>(yWall.fName);
            bc->bname = bname;
            if(bname.find("WALL") != std::string::npos) {
                bc->cname = "DIRICHLET";
                bc->value = Scalar(0);
            } else if(bname.find("interMesh") != std::string::npos) {
                bc->cname = "GHOST";
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
/**
  Find nearest cell
 */
Int Mesh::findNearestCell(const Vector& v) {
    Scalar mindist,dist;
    Int bi = 0;
    mindist = magSq(v - cC[0]);
    for(Int i = 0;i < gBCSfield;i++) {
        dist = magSq(v - cC[i]);
        if(dist < mindist) {
            mindist = dist;
            bi = i;
        }
    }
    return bi;
}
/**
  Find nearest face
 */
Int Mesh::findNearestFace(const Vector& v) {
    Scalar mindist,dist;
    Int bi = 0;
    mindist = magSq(v - fC[0]);
    forEach(fC,i) {
        dist = magSq(v - fC[i]);
        if(dist < mindist) {
            mindist = dist;
            bi = i;
        }
    }
    return bi;
}
/**
  Find nearest cells of probes
 */
void Mesh::getProbeCells(IntVector& probes) {
    forEach(probePoints,j) {
        Vector v = probePoints[j];
        Int index = findNearestCell(v);
        probes.push_back(index);
    }
}
/**
  Find nearest faces of probes
 */
void Mesh::getProbeFaces(IntVector& probes) {
    forEach(probePoints,j) {
        Vector v = probePoints[j];
        Int index = findNearestFace(v);
        probes.push_back(index);
    }
}
/**
  Calculate global courant number
 */
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
/**
  Write all fields
 */
void Mesh::write_fields(Int step) {
    forEachCellField(writeAll(step));
}
/**
  Read all fields
 */
void Mesh::read_fields(Int step) {
    forEachCellField(readAll(step));
}
/**
  Remove all fields
 */
void Mesh::remove_fields() {
    forEachCellField(removeAll());
    forEachFacetField(removeAll());
    forEachVertexField(removeAll());
    BaseField::allFields.clear();
}
/**
  Enroll refine parameters
 */
void Controls::enrollRefine(Util::ParamList& params) {
    params.enroll("direction",&refine_params.dir);
    params.enroll("field",&refine_params.field);
    params.enroll("field_max",&refine_params.field_max);
    params.enroll("field_min",&refine_params.field_min);
    params.enroll("limit",&refine_params.limit);
}
/**
  Enroll domain decomposition parameters
 */
void Controls::enrollDecompose(Util::ParamList& params) {
    params.enroll("n",&decompose_params.n);
    params.enroll("axis",&decompose_params.axis);
    Util::Option* op = new Util::Option(&decompose_params.type, 4, 
            "METIS","XYZ","CELLID","NONE");
    params.enroll("type",op);
}
/**
  Enroll solver control parameters
 */
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
    op = new Option(&write_format,2,"TEXT","BINARY");
    params.enroll("write_format",op);
    op = new Util::BoolOption(&save_average);
    params.enroll("average",op);
    params.enroll("print_time",&print_time);
    params.enroll("npx",&DG::Nop[0]);
    params.enroll("npy",&DG::Nop[1]);
    params.enroll("npz",&DG::Nop[2]);
}
/**
  Create fields
 */
void Prepare::createFields(vector<string>& fields,Int step) {
    BaseField::destroyFields();

    /*for each field*/
    forEach(fields,i) {
        stringstream path;
        Int size;
        path << fields[i] << step;
        std::string str = path.str(); 

        bool exists = false;
        if(System::exists(str+".txt")) {
            exists = true;
            std::ifstream is(str+".txt");
            is >> str >> size;
        } else if(System::exists(str+".bin")) {
            exists = true;
            Util::ifstream_bin is(str+".bin");
            is >> str >> size;
        }

        if(exists) {
            /*fields*/
            BaseField* bf;
            switch(size) {
                case 1 :  bf = new ScalarCellField(fields[i].c_str(),READWRITE,false); break;
                case 3 :  bf = new VectorCellField(fields[i].c_str(),READWRITE,false); break;
                case 6 :  bf = new STensorCellField(fields[i].c_str(),READWRITE,false); break;
                case 9 :  bf = new TensorCellField(fields[i].c_str(),READWRITE,false); break;
            }
            if(bf) bf->fName = bf->fName;
            /*end*/
        }
    }
}
/**
  Read fields
 */
void Prepare::readFields(vector<string>& fields,Int step) {
    Mesh::read_fields(step);
}
/*********************************
 *
 * Adaptive mesh refinement
 *
 *********************************/

/**
  Calculate quantity of interest (QOI)
 */
void Prepare::calcQOI(VectorCellField& qoi) {
    using namespace Mesh;

    BaseField* bf = BaseField::findField(Controls::refine_params.field);
    if(bf) {
        ScalarCellField norm;
        bf->norm(&norm);
        fillBCs(norm);
        applyExplicitBCs(norm,false,false);
        qoi = gradf(norm);
        forEach(qoi,i) {
            Vector& v = qoi[i];
            v[0] = mag(v[0]);
            v[1] = mag(v[1]);
            v[2] = mag(v[2]);
        }
    }
}
/**
  Init AMR refinement threshold
 */
void Prepare::initRefineThreshold() {
    using namespace Mesh;
    using namespace Controls;

    VectorCellField qoi;
    ScalarCellField qoim;
    calcQOI(qoi);
    qoim = mag(qoi);
    /*max and min*/
    Scalar maxq = 0;
    for(Int i = 0;i < gBCS;i++) {
        if(qoim[i] > maxq)
            maxq = qoim[i];
    }
    refine_params.field_max *= (maxq);
    refine_params.field_min *= (maxq);
    std::cout << "----------------------" << std::endl;
    std::cout << " Refinement threshold " << Controls::refine_params.field_max << std::endl;
    std::cout << " Coarsening threshold " << Controls::refine_params.field_min << std::endl;
    std::cout << "----------------------" << std::endl;
}
/**
  Refine grid
 */
void Prepare::refineMesh(Int step,bool init_threshold) {
    using namespace Mesh;
    using namespace Controls;

    std::cout << "Refining grid at step " << step << std::endl;

    /*Specify AMR direction (3D or 2D)*/
    Mesh::amr_direction = refine_params.dir;

    /*Load mesh*/
    LoadMesh(step,false);

    /*create fields*/
    Prepare::createFields(BaseField::fieldNames,step);
    Prepare::readFields(BaseField::fieldNames,step);

    /*init threshold*/
    if(init_threshold)
        Prepare::initRefineThreshold();

    /*find cells to refine/coarsen*/
    IntVector rCells,cCells,rLevel,rDirs;
    {
        /*find quantity of interest*/
        VectorCellField qoi;
        calcQOI(qoi);

        /*get cells to refine*/
        gCells.erase(gCells.begin() + gBCS,gCells.end());
        cCells.assign(gCells.size(),0);
        for(Int i = 0;i < gBCS;i++) {
            Vector q = qoi[i];
            Scalar qm = mag(q);
            RefineParams& rp = refine_params;
            if(qm >= rp.field_max) {
                if(gCells.size() <= rp.limit) {
                    Int dir = 0;
                    if(q[0] >= rp.field_max) dir |= 1;
                    if(q[1] >= rp.field_max) dir |= 2;
                    if(q[2] >= rp.field_max) dir |= 4;

                    if(dir) {
                        Int level = 1;
                        for(;level < 3;level++) {
                            Int factor = (1 << (2 * level));
                            if(qm < rp.field_max * factor)
                                break;
                        }
                        rCells.push_back(i);
                        rLevel.push_back(level);
                        rDirs.push_back(dir);
                    }
                }
            } else if(qm <= rp.field_min) {
                cCells[i] = 1;
            }
        }
    }
    /*Read amr tree*/
    {
        stringstream path;
        int stepn = findLastRefinedGrid(step);
        path << "amrTree" << "_" << stepn << ".bin";
        if(System::exists(path.str())) {
            Util::ifstream_bin is(path.str());
            is >> gAmrTree;
        } else {
            gAmrTree.resize(gCells.size());
            forEach(gCells,i)
                gAmrTree[i].id = i;
        }
    }

    /*refine/coarsen mesh and fields*/
    {
        IntVector refineMap,coarseMap;
        gMesh.refineMesh(rCells,cCells,rLevel,rDirs,refineMap,coarseMap);
        forEachIt(BaseField::allFields, it)
            (*it)->refineField(step,refineMap,coarseMap);
    }
    /*Write amrTree*/
    {
        stringstream path;
        path << "amrTree" << "_" << step << ".bin";
        Util::ofstream_bin os(path.str());
        os << gAmrTree;
    }

    /*Write mesh*/
    {
        stringstream path;
        path << gMeshName << "_" << step << ".bin";
        Util::ofstream_bin os(path.str());
        os << gMesh;
    }

    /*destroy*/
    BaseField::destroyFields();
}

/*********************************
 *
 * Domain Decomposition
 *
 *********************************/
namespace Prepare {

    /**
      Decompose by cell ID
     */
    void decomposeIndex(Int total,IntVector& blockIndex) {
        using namespace Mesh;

        for(Int i = 0;i < gBCS;i++)
            blockIndex[i] = (i / (gBCS / total));
    }
    /**
      Decompose in cartesian directions
     */
    void decomposeXYZ(Int* n,Scalar* nq,IntVector& blockIndex) {
        using namespace Mesh;

        Int i,j,ID;
        Vector maxV(Scalar(-10e30)),minV(Scalar(10e30)),delta;
        Vector axis(nq[0],nq[1],nq[2]);
        Scalar theta = nq[3];
        Vector C;

        /*max and min points*/
        forEach(gVertices,i) {
            C = rotate(gVertices[i],axis,theta);
            for(j = 0;j < 3;j++) {
                if(C[j] > maxV[j]) maxV[j] = C[j];
                if(C[j] < minV[j]) minV[j] = C[j];
            }
        }
        delta = maxV - minV;
        for(j = 0;j < 3;j++) 
            delta[j] /= Scalar(n[j]);

        /*assign block indices to cells*/
        for(i = 0;i < gBCS;i++) {
            C = rotate(gCC[i],axis,theta);
            C = (C - minV) / delta;
            ID = Int(C[0]) * n[1] * n[2] + 
                Int(C[1]) * n[2] + 
                Int(C[2]);
            blockIndex[i] = ID;
        }
    }
    /**
      Decompose using METIS 5.0
     */
    void decomposeMetis(int total,IntVector& blockIndex) {
        using namespace Mesh;

        int ncon = 1;
        int edgeCut = 0;
        int ncells = gBCS;
        std::vector<int> xadj,adjncy;
        std::vector<int> options(METIS_NOPTIONS);

        /*default options*/
        METIS_SetDefaultOptions(&options[0]);

        /*build adjacency*/
        for(Int i = 0;i < gBCS;i++) {
            Cell& c = gCells[i];
            xadj.push_back(adjncy.size());
            forEach(c,j) {
                Int f = c[j];
                if(i == gFOC[f]) {
                    if(gFNC[f] < gBCS)
                        adjncy.push_back(gFNC[f]);
                } else {
                    if(gFOC[f] < gBCS)
                        adjncy.push_back(gFOC[f]);
                }
            }
        }
        xadj.push_back(adjncy.size());

        /*partition*/
        METIS_PartGraphKway (
                &ncells,
                &ncon,
                &xadj[0],
                &adjncy[0],
                NULL,
                NULL,
                NULL,
                &total,
                NULL,
                NULL,
                &options[0],
                &edgeCut,
                (int*)(&blockIndex[0])
                );
    }

}
/**
  Decompose
 */
int Prepare::decomposeMesh(Int step) {
    using namespace Mesh;
    using namespace Constants;

    /*Decompose mesh with master rank*/
    if(MP::host_id == 0) {

        Int total = MP::n_hosts;
        DecomposeParams& dp = Controls::decompose_params;
        vector<string>& fields = BaseField::fieldNames;
        Int i,j,ID,count;

        std::cout << "Decomposing grid at step " << step << std::endl;

        /*Read mesh*/
        LoadMesh(step,false);

        /*Read fields*/
        createFields(fields,step);
        readFields(fields,step);

        /**********************
         * decompose mesh
         **********************/
        std::vector<MeshObject> meshes(total);
        std::map<std::pair<Int,Int>,IntVector> imesh;
        std::vector<IntVector> vLoc(total);
        std::vector<IntVector> fLoc(total);
        std::vector<IntVector> cLoc(total);

        for(i = 0;i < total;i++) {
            vLoc[i].assign(gVertices.size(),0);
            fLoc[i].assign(gFacets.size(),0);
        }

        /*decompose cells*/
        MeshObject *pmesh;
        IntVector *pvLoc,*pfLoc,blockIndex;
        blockIndex.assign(gBCS,0);

        /*choose*/
        if(dp.type == 1) {
            Int tn = dp.n[0] * dp.n[1] * dp.n[2];
            if(total != tn) {
                std::cerr << "Error in XYZ decomposition: use "
                    << tn << " ranks" << std::endl;
                exit(1);
            }
            decomposeXYZ(&dp.n[0],&dp.axis[0],blockIndex);
        } else if(dp.type == 2)
            decomposeIndex(total,blockIndex);
        else
            decomposeMetis(total,blockIndex);

        /***************************************
         * Prepare decomposed grid
         **************************************/

        /*add cells*/
        for(i = 0;i < gBCS;i++) {
            Cell& c = gCells[i];

            /* add cell */
            ID = blockIndex[i];
            pmesh = &meshes[ID];
            pvLoc = &vLoc[ID];
            pfLoc = &fLoc[ID];
            pmesh->mCells.push_back(c);
            cLoc[ID].push_back(i);

            /* mark vertices and facets */
            forEach(c,j) {
                Facet& f = gFacets[c[j]];
                (*pfLoc)[c[j]] = 1;
                forEach(f,k) {
                    (*pvLoc)[f[k]] = 1; 
                }
            }
        }

        /*add vertices & facets*/
        for(ID = 0;ID < total;ID++) {
            pmesh = &meshes[ID];
            pvLoc = &vLoc[ID];
            pfLoc = &fLoc[ID];

            pmesh->mBCS = pmesh->mCells.size();

            count = 0;
            forEach(gVertices,i) {
                if((*pvLoc)[i]) {
                    pmesh->mVertices.push_back(gVertices[i]);
                    (*pvLoc)[i] = count++;
                } else
                    (*pvLoc)[i] = Constants::MAX_INT;
            }

            count = 0;
            forEach(gFacets,i) {
                if((*pfLoc)[i]) {
                    pmesh->mFacets.push_back(gFacets[i]);
                    (*pfLoc)[i] = count++;
                } else
                    (*pfLoc)[i] = Constants::MAX_INT;
            }
        }
        /*adjust IDs*/
        for(ID = 0;ID < total;ID++) {
            pmesh = &meshes[ID];
            pvLoc = &vLoc[ID];
            pfLoc = &fLoc[ID];

            forEach(pmesh->mFacets,i) {
                Facet& f = pmesh->mFacets[i];
                forEach(f,j)
                    f[j] = (*pvLoc)[f[j]];
            }

            forEach(pmesh->mCells,i) {
                Cell& c = pmesh->mCells[i];
                forEach(c,j)
                    c[j] = (*pfLoc)[c[j]];
            }
        }
        /*inter-processor boundary faces*/
        Int co,cn;
        forEach(gFacets,i) {
            if(gFNC[i] < gBCS) {
                co = blockIndex[gFOC[i]];
                cn = blockIndex[gFNC[i]];
                if(co != cn) {
                    imesh[std::make_pair(co,cn)].push_back(fLoc[co][i]);
                    imesh[std::make_pair(cn,co)].push_back(fLoc[cn][i]);
                }
            }
        }
        /***************************
         * Boundaries
         ***************************/
        for(ID = 0;ID < total;ID++) {
            pmesh = &meshes[ID];
            pfLoc = &fLoc[ID];

            /*physical boundary*/
            forEachIt(gMesh.mBoundaries,it) {
                IntVector b;    
                Int f;
                forEach(it->second,j) {
                    f = (*pfLoc)[it->second[j]];
                    if(f != Constants::MAX_INT)
                        b.push_back(f);
                }
                if(b.size())
                    pmesh->mBoundaries[it->first] = b;
            }

            /*inter-processor boundaries*/
            for(j = 0;j < total;j++) {
                std::pair<Int,Int> key = std::make_pair(ID,j);
                if(imesh.find(key) != imesh.end()) {
                    stringstream path;
                    path << "interMesh_" << ID << "_" << j;
                    string str = path.str();
                    pmesh->mBoundaries[str.c_str()] = imesh[key];
                }
            }

        }
        /**************************************
         * Write grid
         **************************************/

        /* write master grid after reordering cells*/
        {
            Cells reordered;
            for(ID = 0;ID < total;ID++) {
                forEach(cLoc[ID],i) {
                    Cell& c = gCells[cLoc[ID][i]];
                    reordered.push_back(c);
                }
            }
            for(i = 0;i < gBCS;i++)
                gCells[i] = reordered[i];

            /*write mesh*/
            stringstream path;
            path << gMeshName << "_" << step;
            string str = path.str();

            if(System::exists(str + ".txt")) {
                ofstream os(str + ".txt");
                os << gMesh;
            } else if(System::exists(str + ".bin")) {
                Util::ofstream_bin os(str + ".bin");
                os << gMesh;
            }
        }

        /* write slave grids*/
        for(ID = 0;ID < total;ID++) {
            pmesh = &meshes[ID];

            /*create directory*/
            stringstream path;
            path << gMeshName << ID;
            System::mkdir(path.str());

            /*mesh*/
            stringstream path1;
            path1 << path.str() << "/" << gMeshName << "_" << step;
            string str = path1.str();

            if(Controls::write_format == Controls::TEXT) {
                ofstream os(str + ".txt");
                os << *pmesh << "\n";
            } else {
                Util::ofstream_bin os(str + ".bin");
                os << *pmesh << "\n";
            }
        }

        /***************************
         * Expand cLoc
         ***************************/
        for(ID = 0;ID < total;ID++) {
            const Int block = DG::NP;
            IntVector& cF = cLoc[ID];
            cF.resize(cF.size() * block);
            for(int i = cF.size() - 1;i >= 0;i -= block) {
                Int ii = i / block;
                Int C = cF[ii] * block;
                for(Int j = 0; j < block;j++)
                    cF[i - j] = C + block - 1 - j; 
            }
        }
        /***************************
         * Write fields
         ***************************/

        /*slave fields*/
        for(ID = 0;ID < total;ID++) {
            pmesh = &meshes[ID];

            forEach(fields,i) {
                BaseField* pf = BaseField::findField(fields[i]);
                if(!pf) continue;

                stringstream path;
                path << gMeshName << ID << "/" << fields[i] << step;
                string str = path.str();

                if(Controls::write_format == Controls::TEXT) {
                    std::ofstream os(str + ".txt");
                    pf->write(os, &cLoc[ID]);
                    /*inter mesh boundaries*/
                    forEachIt(pmesh->mBoundaries,it) {
                        if(it->first.find("interMesh") != string::npos)
                           os << it->first << " " << "{\n\ttype GHOST\n}\n";
                    }
                } else {
                    Util::ofstream_bin os(str + ".bin");
                    pf->write(os, &cLoc[ID]);
                    /*inter mesh boundaries*/
                    forEachIt(pmesh->mBoundaries,it) {
                        if(it->first.find("interMesh") != string::npos)
                           os << it->first << " " << "{\n\ttype GHOST\n}\n";
                    }
                }

            }
        }

        /*master fields*/
        IntVector cLocAll(gBCS * DG::NP);
        count = 0;
        for(ID = 0;ID < total;ID++) {
            forEach(cLoc[ID], i)
                cLocAll[count++] = cLoc[ID][i];
        }
        forEach(fields,i) {
            BaseField* pf = BaseField::findField(fields[i]);
            if(!pf) continue;

            stringstream path;
            path << fields[i] << step;
            string str = path.str();

            if(System::exists(str + ".txt")) {
                ofstream os(str + ".txt");
                pf->write(os, &cLocAll);
            } else if(System::exists(str + ".bin")) {
                Util::ofstream_bin os(str + ".bin");
                pf->write(os, &cLocAll);
            }
        }

        /*destroy*/ 
        BaseField::destroyFields();
    }
    MP::barrier();

    /*change directory*/
    stringstream s;
    s << Mesh::gMeshName << MP::host_id;
    System::cd(s.str());
    Mesh::LoadMesh(step);  

    return 0;
}
/**
  Reverse decomposition
 */
int Prepare::mergeFields(Int step) {
    using namespace Mesh;
    using namespace Controls;

    vector<string>& fields = BaseField::fieldNames;

    /*indexes*/
    Int total = MP::n_hosts;

    std::cout << "Merging fields at step " << step << std::endl;

    /*Read mesh*/
    int stepm = findLastRefinedGrid(step);
    LoadMesh(stepm,false);

    /*Read fields*/
    createFields(fields,stepm);
    readFields(fields,stepm);

    /*read and merge fields*/
    forEach(fields,i) {
        BaseField* pf = BaseField::findField(fields[i]);
        if(!pf) continue;

        Int offset = 0;
        for(Int ID = 0;ID < total;ID++) {
            stringstream fpath;
            fpath << gMeshName << ID << "/" << fields[i] << step;
            string str = fpath.str();

            Int nread;
            if(System::exists(str + ".txt")) {
                std::ifstream is(str + ".txt");
                nread = pf->readInternal(is,offset);
            } else if(System::exists(str + ".bin")) {
                Util::ifstream_bin is(str + ".bin");
                nread = pf->readInternal(is,offset);
            } else {
                continue;
            }
            offset += nread;
        }
    }
    write_fields(step);

    /*destroy*/ 
    BaseField::destroyFields();

    return 0;
}

