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
	
	Util::ParamList params("general");
	params.enroll("mesh",&Mesh::gMeshName);
	params.enroll("fields",&fields);
	params.enroll("decompose",&n);
	params.enroll("axis",&axis);
	params.enroll("probe",&Mesh::probePoints);
	Util::Option* op = new Util::Option(&decomptype, 4, 
			"XYZ","INDEX","METIS","NONE");
	params.enroll("decomptype",op);
	params.read(input); 

	/*Mesh*/
	if(mp.n_hosts > 1) {
		stringstream s;
		s << Mesh::gMeshName << mp.host_id;
		if(!System::cd(s.str()))
			return 1;
	}
	Mesh::readMesh();
	Mesh::initGeomMeshFields(work != 0);
	cout << "fields " << fields << endl;
	atexit(Util::cleanup);

	/*do work*/
	if(work == 1) {
		Prepare::merge(Mesh::gMesh,&n[0],fields,Mesh::gMeshName,start_index);
	} else if(work == 2) {
		Prepare::convertVTK(Mesh::gMesh,fields,start_index);
	} else if(work == 3) {
		Prepare::probe(Mesh::gMesh,fields,start_index);
	} else{
		Prepare::decompose(Mesh::gMesh,&n[0],&axis[0],decomptype);
		Prepare::decomposeFields(fields,Mesh::gMeshName,start_index);
	} 
	return 0;
}
