#include <iostream>
#include <string>
#include "ModelClass.h"

using namespace std;

int main(int argc, char* argv[]) {
// ------------------------------------
//int main() {
//	int argc = 2;
//	string argv[2];
//	argv[1] = "C:/Users/evans/ASenDHome/ASenD3D/examples/testCases/shellBeam/distributedLoading/job.yaml";
// --------------------------------------

	if(argc < 2) {
		cout << "Error: no job script file name provided." << endl;
		return 0;
	} else {
		string job_file = argv[1];
		cout << "Beginning Job: " << job_file << endl;
		Model job_model;
		job_model.read_job(job_file);
		job_model.execute_job();
		return 0;
	}
}