#include "Mesher.h"
#include <string>
#include <iostream>

using namespace std;

//int main(int argc, char* argv[]) {
//--------------------------------
int main() {
	int argc = 3;
	char myChAr[20];
	string argv[3];
	argv[1] = "C:/Users/evans/ASenDHome/surfaceMesh.yaml";
	argv[2] = "C:/Users/evans/ASenDHome/tetMesh.yaml";
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
	mesher.generateMesh();
	mesher.distributeNodes();
	mesher.writeOutput(outputFile);
	return 0;
}