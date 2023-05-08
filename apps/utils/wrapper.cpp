#include "field.h"
#include "mp.h"
#include "system.h"
#include "wrapper.h"

namespace Solver {

    std::ifstream input;
    static std::string sname;

    void Initialize(int argc, char* argv[]) {
        /*cleanup*/
        atexit(MP::cleanup);
        /*message passing object*/
        MP::printOn = (MP::host_id == 0);
        if(!strcmp(argv[1],"-h")) {
            std::cout << "Usage:\n"
                << "  " << argv[0] << " <inputfile>\n"
                << "Options:\n"
                << "  -h          --  Display this message\n\n";
            exit(0);
        } 
        input.open(argv[1]);

        /*General options*/
        {
            Util::ParamList params("general");
            params.enroll("solver", &sname);
            params.enroll("mesh", &Mesh::gMeshName);
            Mesh::enroll(params);
            params.read(input);
            std::string path(argv[0]);
            if(path.find(sname) == std::string::npos) {
                std::cout << "Incorrect solver " << path
                          << " used instead of " << sname
                          << "." << std::endl;
                exit(1);
            }
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
    }

    void Finalize() {
        input.close();
#ifdef _DEBUG
        /*print memory usage*/
        std::cout << "====================================" << std::endl;
        std::cout << "Memory Usage:" << std::endl;
        forEachCellField(printUsage());
        forEachFacetField(printUsage());
        forEachVertexField(printUsage());
        std::cout << "====================================" << std::endl;
#endif  
    }
}
