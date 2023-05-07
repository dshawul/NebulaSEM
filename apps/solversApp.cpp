#include "field.h"
#include "mp.h"
#include "system.h"
#include "solve.h"
#include "iteration.h"

#include "walldist.h"
#include "potential.h"
#include "convection.h"
#include "diffusion.h"
#include "transport.h"
#include "wave.h"
#include "piso.h"
#include "euler.h"
#include "hydro_balance.h"

using namespace std;

/**
  \verbatim
  Main application entry point for different solvers.
  \endverbatim
 */
int main(int argc, char* argv[]) {

    /*message passing object*/
    MP mp(argc, argv);
    MP::printOn = (MP::host_id == 0);
    if(!strcmp(argv[1],"-h")) {
        std::cout << "Usage:\n"
            << "  ./solvers <inputfile>\n"
            << "Options:\n"
            << "  -h          --  Display this message\n\n";
        return 0;
    } 
    ifstream input(argv[1]);

    /*General options*/
    string sname;
    {
        Util::ParamList params("general");
        params.enroll("solver", &sname);
        params.enroll("mesh", &Mesh::gMeshName);
        Mesh::enroll(params);
        params.read(input);
    }
    /*AMR options*/
    {
        Util::ParamList params("refinement");
        Controls::enrollRefine(params);
        params.read(input);
    }
    /*Decompose options*/
    {
        Util::ParamList params("decomposition");
        Controls::enrollDecompose(params);
        params.read(input); 
    }
    /*Fields*/
    {
        Util::ParamList params("prepare");
        params.enroll("fields",&BaseField::fieldNames);
        params.read(input);
    }
    /*cleanup*/
    atexit(MP::cleanup);

    /*call solver*/
    if (!Util::compare(sname, "piso")) {
        piso(input);
    } else if (!Util::compare(sname, "euler")) {
        euler(input);
    } else if (!Util::compare(sname, "diffusion")) {
        diffusion(input);
    } else if (!Util::compare(sname, "convection")) {
        convection(input);
    } else if (!Util::compare(sname, "transport")) {
        transport(input);
    } else if (!Util::compare(sname, "potential")) {
        potential(input);
    } else if (!Util::compare(sname, "hydro_balance")) {
        hydro_balance(input);
    } else if (!Util::compare(sname, "walldist")) {
        walldist(input);
    } else if (!Util::compare(sname, "wave")) {
        wave(input);
    }

#ifdef _DEBUG
    /*print memory usage*/
    std::cout << "====================================" << std::endl;
    std::cout << "Memory Usage:" << std::endl;
    forEachCellField(printUsage());
    forEachFacetField(printUsage());
    forEachVertexField(printUsage());
    std::cout << "====================================" << std::endl;
#endif  

    return 0;
}
