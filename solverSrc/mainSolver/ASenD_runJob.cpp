#include <iostream>
#include <string>
#include "ModelClass.h"

using namespace std;

int main(int argc, char* argv[]) {
// ------------------------------------
//int main() {
//	int argc = 2;
//	string argv[2];
//	argv[1] = "C:/Users/evans/ASenDHome/ASenD3D/examples/dynamicAnalysis/shellDiskImpact/job.yaml";
// --------------------------------------

	if(argc < 2) {
		cout << "Error: no job script file name provided." << endl;
		return 0;
	} else {
		string jobFile = argv[1];
		cout << "Beginning Job: " << jobFile << endl;
		Model jobModel;
		jobModel.readJob(jobFile);
		jobModel.executeJob();
		return 0;
	}
}