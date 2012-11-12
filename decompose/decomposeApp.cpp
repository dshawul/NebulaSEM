#include "mesh.h"
#include "decompose.h"

using namespace std;

/*decompose application*/
int main(int argc,char* argv[]) {
	/*input stream*/
	ifstream input(argv[1]);

    /*read mesh & fields to decompose*/
	input >> Mesh::gMeshName;
	Mesh::readMesh();
	IntVector n;
	input >> n;
	vector<string> fields;
	input >> fields;
	cout << "fields " << fields << endl;

	/*do XYZ decomposition*/
	Decompose::decomposeXYZ(Mesh::gMesh,&n[0]);
	Decompose::decomposeFields(fields,Mesh::gMeshName);
	return 0;
}
