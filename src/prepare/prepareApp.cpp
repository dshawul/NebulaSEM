#include "mesh.h"
#include "prepare.h"
#include "system.h"

using namespace std;

/**
\verbatim
Post/Pre processing jobs such as
  Domain decomposition
  Merging results of decomposed domains
  Converting data to VTK format
  Probing result at specified locations
\endverbatim
*/
int main(int argc,char* argv[]) {
    /*message passing object*/
    MP mp(argc,argv);
    MP::printOn = (MP::host_id == 0);
    /*cmd line*/
    int work = 0;
    Int start_index = 0;
    for(int i = 1;i < argc;i++) {
        if(!strcmp(argv[i],"-merge")) {
            work = 1;
        } else if(!strcmp(argv[i],"-vtk")) {
            work = 2;
        } else if(!strcmp(argv[i],"-probe")) {
            work = 3;
        } else if(!strcmp(argv[i],"-refine")) {
            work = 4;
        } else if(!strcmp(argv[i],"-poly")) {
            Vtk::write_polyhedral = true;
        } else if(!strcmp(argv[i],"-start")) {
            i++;
            start_index = atoi(argv[i]);
        } else if(!strcmp(argv[i],"-h")) {
            std::cout << "Usage:\n"
                      << "  ./prepare <inputfile> <Options>\n"
                      << "Options:\n"
                      << "  -merge      --  Merge results of decomposed domain\n"
                      << "  -vtk        --  Convert data to VTK format\n"
                      << "  -probe      --  Probe result at specified locations\n"
                      << "  -refine     --  Refine mesh\n"
                      << "  -poly       --  Write VTK in polyhedral format\n"
                      << "  -start <i>  --  Start at time step <i>\n"
                      << "  -h          --  Display this message\n\n";
            return 0;
        } 
    }
    
    /*open job specification file*/
    ifstream input(argv[1]);

    /*read mesh & fields*/
    vector<string>& fields = BaseField::fieldNames;
    
    /*Options*/
    {
        Util::ParamList params("general");
        params.enroll("mesh",&Mesh::gMeshName);
        params.enroll("probe",&Mesh::probePoints);
        params.enroll("npx",&DG::Nop[0]);
        params.enroll("npy",&DG::Nop[1]);
        params.enroll("npz",&DG::Nop[2]);
        params.enroll("gravity",&Controls::gravity);
        params.read(input); 
    }
    {
        Util::ParamList params("prepare");
        params.enroll("fields",&fields);
        params.read(input);
    }
    {
        Util::ParamList params("decomposition");
        Controls::enrollDecompose(params);
        params.read(input); 
    }
    {
        Util::ParamList params("refinement");
        Controls::enrollRefine(params);
        params.read(input);
    }

    /*switch directory*/
    if(mp.n_hosts > 1) {
        stringstream s;
        s << Mesh::gMeshName << mp.host_id;
        if(!System::cd(s.str()))
            return 1;
    }
    atexit(MP::cleanup);

    /*do work*/
    if(MP::printOn)
        cout << "--------------------------------------------\n";
    if(work == 1) {
        cout << "Merging decomposed domain.\n";
        Prepare::mergeFields(start_index);
    } else if(work == 2) {
        cout << "Converting result to VTK format.\n";
        Prepare::convertVTK(fields,start_index);
    } else if(work == 3) {
        cout << "Probing result at specified locations.\n";
        Prepare::probe(fields,start_index);
    } else if(work == 4) {
        cout << "Refining grid.\n";
        Prepare::refineMesh(start_index);
    } else {
        cout << "Decomposing domain.\n";
        Prepare::decomposeMesh(start_index);
    } 
    return 0;
}
