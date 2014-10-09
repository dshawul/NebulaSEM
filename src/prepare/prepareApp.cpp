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
	
	/*cmd line*/
	int work = 0;
	int refine_type = 0,refine_shape = 0;
	Vector refine_dir(0.0);
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
    Int decomptype = 2;
	vector<string> fields;
	vector<Int> n;
	vector<Scalar> axis(4);
	axis[0] = 1;
	
	/*Options*/
	{
		Util::ParamList params("general");
		params.enroll("mesh",&Mesh::gMeshName);
		params.enroll("fields",&fields);
		params.enroll("probe",&Mesh::probePoints);
		params.read(input); 
	}
	{
		Util::ParamList params("decomposition");
		params.enroll("n",&n);
		params.enroll("axis",&axis);
		Util::Option* op = new Util::Option(&decomptype, 4, 
				"XYZ","INDEX","METIS","NONE");
		params.enroll("decomptype",op);
		params.read(input); 
	}
	{
		Util::ParamList params("refinement");
		Util::Option* op;
		op = new Util::Option(&refine_type, 2,"CELL","FACE");
		params.enroll("type",op);
		op = new Util::Option(&refine_shape, 2,"QUAD","TRI");
		params.enroll("shape",op);
		params.enroll("direction",&refine_dir);
		params.read(input); 
	}
	/*Mesh*/
	if(mp.n_hosts > 1) {
		stringstream s;
		s << Mesh::gMeshName << mp.host_id;
		if(!System::cd(s.str()))
			return 1;
	}
	Mesh::readMesh();
	Mesh::initGeomMeshFields(false);
	cout << "fields " << fields << endl;
	atexit(Util::cleanup);

	/*do work*/
	if(work == 1) {
		Prepare::merge(&n[0],fields,Mesh::gMeshName,start_index);
	} else if(work == 2) {
		Prepare::convertVTK(fields,start_index);
	} else if(work == 3) {
		Prepare::probe(fields,start_index);
	} else if(work == 4) {
		IntVector rCells;
		for(Int i = 0;i < Mesh::gBCellsStart;i++)
			rCells.push_back(i);
		Prepare::refineMesh(rCells,refine_type,refine_shape,refine_dir);
	} else {
		Prepare::decompose(&n[0],&axis[0],decomptype);
		Prepare::decomposeFields(fields,Mesh::gMeshName,start_index);
	} 
	return 0;
}
