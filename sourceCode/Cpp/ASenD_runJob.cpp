#include <iostream>
#include <string>
#include "ModelClass.h"

using namespace std;

int main(int argc, char* argv[]) {
	
	if(argc < 2) {
		cout << "Error: no job script file name provided." << endl;
		return 0;
	} else {
		string jobFile = argv[1];
		cout << "Beginning Job: " << jobFile << endl;
		Model jobModel;
		jobModel.readJob(jobFile);
		jobModel.executeJob();
		jobModel.destroy();
		return 0;
	}
}