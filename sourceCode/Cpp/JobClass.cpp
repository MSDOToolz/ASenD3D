#include "JobClass.h"
#include <string>

using namespace std;

JobCommand::JobCommand() {
	cmdString = "";
	fileName = "";
	
	// for solve
	dynamic = false;
	elastic = true;
	loadRampSteps = 1;
	newmarkBeta = 0.25;
	newmarkGamma = 0.5;
	nonlinearGeom = false;
	saveSolnHist = false;
	simPeriod = 1.0;
	solverBandwidth = 100;
	solverBlockDim = 100;
	solverMethod = "direct";
	staticLoadTime = 0.0;
	thermal = false;
	timeStep = 1.0;
	
	// modal
	type = "buckling";
	numModes = 10;
	
	// setSolnToMode
	solnField = "displacement";
	mode = 1;
	maxAmplitude = 1.0;
	
	// writeNodeResults
	nodeSet = "all";
	
	// writeElementResults
	elementSet = "all";
	
	// writeModalResults
	writeModes = true;
	
	// writeObjective
	writeGradient = true;	
	
	next = NULL;
	return;
}

Job::Job() {
	firstCmd = NULL;
	lastCmd = NULL;
	length = 0;
}

void Job::addCommand(JobCommand *newCmd) {
	if(!firstCmd) {
		firstCmd = newCmd;
		lastCmd = newCmd;
	} else {
		lastCmd->next = newCmd;
		lastCmd = newCmd;
	}
	length++;
}

int Job::getLength() {
	return length;
}

JobCommand* Job::getFirst() {
	return firstCmd;
}