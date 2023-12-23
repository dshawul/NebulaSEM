#include "prepare.h"

using namespace std;
using namespace Mesh;

/**
  Convert to VTK format
 */
void Prepare::convertVTK(std::vector<std::string>& fields, Int start_index, Int stop_index) {
    for(Int step = start_index; step < stop_index; step++) {
        if(LoadMesh(step,false))
            createFields(fields,step);
        readFields(fields,step);
        /*write vtk*/
        Vtk::write_vtk(step);
    }
}
/**
  Probe values at specified locations
 */
void Prepare::probe(std::vector<std::string>& fields, Int start_index, Int stop_index) {
    /*probe points*/
    IntVector probes;
    ofstream of("probes");

    /*probe at each time step*/
    for(Int step = start_index; step < stop_index; step++) {
        if(LoadMesh(step,true)) {
            createFields(fields,step);
            probes.clear();
            getProbeFaces(probes);
        }
        readFields(fields,step);

        /*Interpolate*/
        forEachCellField(interpolateVertexAll());

        /*write probes*/
#define ADD(v,value,weight) {                                   \
    dist = magSq((v) - probeP);                                 \
    dist = weight / (dist + 1.0f);                              \
    sum += (value) * dist;                                      \
    sumd += dist;                                               \
}
#define SUM(X) {                                                \
    Cell& c = gCells[X];                                        \
    forEach(c,m) {                                              \
        Facet& f = gFacets[c[m]];                               \
        forEach(f,j) {                                          \
            ADD(gVertices[f[j]],(*it)[f[j]],1.0);               \
        }                                                       \
    }                                                           \
}
#define WRITE(T) {                                              \
    auto it1 = MeshField<T,CELL>::fields_.begin();              \
    auto vff = MeshField<T,CELL>::vf_fields_;                   \
    forEachIt(*vff,it) {                                        \
        T sum(0.0);                                             \
        Scalar sumd(0.0);                                       \
        ADD(cC[c1],(*(*it1))[c1],2.0);                          \
        ADD(cC[c2],(*(*it1))[c2],2.0);                          \
        SUM(sc);                                                \
        of <<  (sum/sumd) << " ";                               \
        it1++;                                                  \
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

            of << "\n";
        }
#undef WRITE
#undef SUM
#undef ADD
    }
}


