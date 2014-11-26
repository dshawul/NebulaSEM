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
    RefineParams ref_params;
    Int decomptype = 3;
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
		params.enroll("type",op);
		params.read(input); 
	}
	{
		Util::ParamList params("refinement");
		Util::Option* op;
		op = new Util::Option(&ref_params.shape, 2,"QUAD","TRI");
		params.enroll("shape",op);
		params.enroll("direction",&ref_params.dir);
		params.enroll("field",&ref_params.field);
		params.enroll("field_max",&ref_params.field_max);
		params.enroll("field_min",&ref_params.field_min);
		params.enroll("limit",&ref_params.limit);
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
	if(work == 1) {
		cout << "Merging decomposed domain.\n";
		Prepare::merge(&n[0],fields,start_index);
	} else if(work == 2) {
		cout << "Converting result to VTK format.\n";
		Prepare::convertVTK(fields,start_index);
	} else if(work == 3) {
		cout << "Probing result at specified locations.\n";
		Prepare::probe(fields,start_index);
	} else if(work == 4) {
		cout << "Refining grid.\n";
		Prepare::refineMesh(fields,ref_params,start_index);
	} else {
		cout << "Decomposing domain.\n";
		Prepare::decompose(fields,&n[0],&axis[0],decomptype,start_index);
	} 
	return 0;
}
