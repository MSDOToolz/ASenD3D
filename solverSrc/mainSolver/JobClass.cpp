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
	rayDampCK = 0.0;
	rayDampCM = 0.0;
	nonlinearGeom = false;
	saveSolnHist = false;
	lumpMass = false;
	fullReform = 1;
	simPeriod = 1.0;
	solverBandwidth = 2000000000;
	solverBlockDim = 2000000000;
	solverMethod = "direct";
	maxIt = 0;
	convTol = 1.0e-12;
	thermal = false;
	timeStep = 1.0;
	
	// modal
	type = "buckling";
	numModes = 10;
	tgtEval = 0.0;
	
	// setSolnToMode
	solnField = "displacement";
	mode = 1;
	maxAmplitude = 1.0;
	
	// writeNodeResults
	nodeSet = "all";
	timeStepTag = "";
	
	// writeElementResults
	elementSet = "all";
	position = "centroid";
	
	// writeModalResults
	writeModes = true;
	
	// writeObjective
	writeGradient = true;	
	return;
}

