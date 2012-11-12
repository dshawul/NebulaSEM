#include "vtk.h"

using namespace std;

/*vtk application*/
int main(int argc,char* argv[]) {
	/*input stream*/
	ifstream input(argv[1]);

    /*read mesh*/
	input >> Mesh::gMeshName;
	Mesh::readMesh();
	vector<string> fields;
	input >> fields;
	cout << "fields " << fields << endl;

	return 0;
}
