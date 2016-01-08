#include "prepare.h"

using namespace std;
using namespace Mesh;

/**
Convert to VTK format
*/
int Prepare::convertVTK(vector<string>& fields,Int start_index) {
    for(Int step = start_index;;step++) {
        if(LoadMesh(step,(step == start_index),true))
            createFields(fields,step);
        if(!readFields(fields,step))
            break;
        /*write vtk*/
        Vtk::write_vtk(step);
    }

    return 0;
}
/**
Probe values at specified locations
*/
int Prepare::probe(vector<string>& fields,Int start_index) {
    /*probe points*/
    IntVector probes;
    ofstream of("probes");

    /*probe at each time step*/
    for(Int step = start_index;;step++) {
        if(LoadMesh(step,(step == start_index),true)) {
            createFields(fields,step);
            probes.clear();
            getProbeFaces(probes);
        }
        if(!readFields(fields,step))
            break;

        /*Interpolate*/
        forEachCellField(interpolateVertexAll());
        
        /*write probes*/
#define ADD(v,value,weight) {                                       \
        dist = magSq((v) - probeP);                                 \
        dist = weight / (dist + 1.0f);                              \
        sum += (value) * dist;                                      \
        sumd += dist;                                               \
}
#define SUM(X) {                                                    \
        Cell& c = gCells[X];                                        \
        forEach(c,m) {                                              \
            Facet& f = gFacets[c[m]];                               \
            forEach(f,j) {                                          \
                ADD(gVertices[f[j]],(*it)[f[j]],1.0);               \
            }                                                       \
        }                                                           \
}
#define WRITE(T) {                                                  \
        std::list<MeshField<T,CELL>*>::iterator it1 =               \
            MeshField<T,CELL>::fields_.begin();                     \
        for(MeshField<T,CELL>::vertexFieldsType::iterator it =      \
            (MeshField<T,CELL>::vf_fields_)->begin(); it !=         \
            (MeshField<T,CELL>::vf_fields_)->end(); ++it,++it1) {   \
            T sum(0.0);                                             \
            Scalar sumd(0.0);                                       \
            ADD(cC[c1],(*(*it1))[c1],2.0);                          \
            ADD(cC[c2],(*(*it1))[c2],2.0);                          \
            SUM(sc);                                                \
            of <<  (sum/sumd) << " ";                               \
        }                                                           \
}
        forEach(probes,i) {
            Int fi = probes[i];
            Int c1 = FO[fi];
            Int c2 = FN[fi];
            Vector probeP = probePoints[i];
            Scalar dir = ((fC[fi] - probeP) & fN[fi]),dist;
            Int sc = (dir >= 0) ? c1 : c2;

            of << step << " " << i << " " << probePoints[i] << " ";

            WRITE(Scalar);
            WRITE(Vector);
            WRITE(STensor);
            WRITE(Tensor);

            of << endl;
        }
#undef WRITE
#undef SUM
#undef ADD
    }

    return 0;
}


