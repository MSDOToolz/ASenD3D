#include "Mesher.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
//--------------------------------
//int main() {
//	int argc = 3;
//	string argv[3];
//	argv[1] = "C:/Users/evans/ASenDHome/ASenD3D/examples/testCases/tetBeam/unstructInp.yaml";
//	argv[2] = "C:/Users/evans/ASenDHome/ASenD3D/examples/testCases/tetBeam/tetMesh.yaml";
//-------------------------------------
	Mesher mesher;
	string inputFile;
	string outputFile;
	if (argc > 2) {
		inputFile = argv[1];
		outputFile = argv[2];
	}
	else if (argc > 1) {
		inputFile = argv[1];
		outputFile = "meshOut.yaml";
	}
	else {
		cout << "Error: no input file provided for unstructured tet mesh." << endl;
		return 1;
	}
	mesher.readInput(inputFile);
	mesher.prep();
	cout << "generating mesh" << endl;
	bool resolved = mesher.generateMesh();
	cout << "finished mesh" << endl;
	if (resolved) {
		mesher.distributeNodes();
		cout << "distributing nodes" << endl;
	}
	else {
		cout << "Warning: mesh reached maximum number of elements unresolved." << endl;
	}
	mesher.writeOutput(outputFile);
	return 0;
}