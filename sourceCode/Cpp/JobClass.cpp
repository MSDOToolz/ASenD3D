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
	solverBlockDim = 2000000000;
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
	timeStepTag = "";
	
	// writeElementResults
	elementSet = "all";
	
	// writeModalResults
	writeModes = true;
	
	// writeObjective
	writeGradient = true;	
	
	next = nullptr;
	return;
}

void JobCommand::destroy() {
	fields.destroy();
	timeSteps.destroy();
	properties.destroy();
	objInclude.destroy();
	return;
}

Job::Job() {
	firstCmd = nullptr;
	lastCmd = nullptr;
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

void Job::destroy() {
	JobCommand* thisCmd = firstCmd;
	JobCommand* nextCmd;
	while (thisCmd) {
		nextCmd = thisCmd->next;
		thisCmd->destroy();
		delete thisCmd;
		thisCmd = nextCmd;
	}
	return;
}