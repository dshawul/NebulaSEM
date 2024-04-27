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
    Int start_index = 0, stop_index = 0;
    for(int i = 1;i < argc;i++) {
        if(!strcmp(argv[i],"-merge")) {
            work = 1;
        } else if(!strcmp(argv[i],"-vtk")) {
            work = 2;
        } else if(!strcmp(argv[i],"-probe")) {
            work = 3;
        } else if(!strcmp(argv[i],"-refine")) {
            work = 4;
        } else if(!strcmp(argv[i],"-coords")) {
            work = 5;
        } else if(!strcmp(argv[i],"-start")) {
            i++;
            start_index = atoi(argv[i]);
            stop_index = start_index + 1;
        } else if(!strcmp(argv[i],"-stop")) {
            i++;
            stop_index = atoi(argv[i]);
        } else if(!strcmp(argv[i],"-h")) {
            std::cout << "Usage:\n"
                << "  ./prepare <inputfile> <Options>\n"
                << "Options:\n"
                << "  -merge      --  Merge results of decomposed domain\n"
                << "  -vtk        --  Convert data to VTK format\n"
                << "  -probe      --  Probe result at specified locations\n"
                << "  -refine     --  Refine mesh\n"
                << "  -coords     --  Write coordinates to a file\n"
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
        Util::Option* op;
        params.enroll("gravity", &Controls::gravity);
        op = new Util::BoolOption(&Mesh::is_spherical);
        params.enroll("is_spherical",op);
        params.enroll("sphere_radius", &Mesh::sphere_radius);
        op = new Util::Option(&Controls::write_format,{"TEXT","BINARY"});
        params.enroll("write_format",op);
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
    {
        Util::ParamList params("vtk");
        Util::Option* op;
        op = new Util::BoolOption(&Vtk::write_polyhedral);
        params.enroll("write_polyhedral",op);
        op = new Util::BoolOption(&Vtk::write_cell_value);
        params.enroll("write_cell_value",op);
        params.read(input);
    }

    /*switch directory*/
    if(work != 0 && mp.n_hosts > 1) {
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
        if(MP::host_id == 0 && MP::n_hosts > 1) {
            cout << "Merging decomposed domain.\n";
            System::cd(MP::workingDir);
            Prepare::mergeFields(start_index);
        }
    } else if(work == 2) {
        cout << "Converting result to VTK format.\n";
        Prepare::convertVTK(fields,start_index,stop_index);
    } else if(work == 3) {
        cout << "Probing result at specified locations.\n";
        Prepare::probe(fields,start_index,stop_index);
    } else if(work == 4) {
        if(MP::host_id == 0) {
            cout << "Refining grid.\n";
            Prepare::refineMesh(start_index);
        }
    } else if(work == 5) {
        cout << "Writing coordinates.\n";
        Prepare::writeCoords(start_index);
    } else {
        System::cd(MP::workingDir);
        if(MP::n_hosts > 1)
            Prepare::decomposeMesh(start_index);
    } 
    MP::barrier();
    return 0;
}
