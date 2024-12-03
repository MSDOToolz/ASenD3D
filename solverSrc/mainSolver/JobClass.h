#ifndef JOBCLASS
#define JOBCLASS
#include "ListEntClass.h"
#include <vector>
#include <string>

class JobCommand {
	public:
	    std::string cmdString;
		std::string fileName;
		
		// for solve
		bool dynamic;
		bool elastic;
		int loadRampSteps;
		double newmarkBeta;
		double newmarkGamma;
		double rayDampCK;
		double rayDampCM;
		bool nonlinearGeom;
		bool saveSolnHist;
		double simPeriod;
		bool lumpMass;
		int fullReform;
		int solverBandwidth;
		int solverBlockDim;
		std::string solverMethod;
		int maxIt;
		double convTol;
		std::list<double> staticLoadTime;
		bool thermal;
		double timeStep;
		
		// modal
		std::string type;
		int numModes;
		double tgtEval;
		
		// setSolnToMode
		std::string solnField;
		int mode;
		double maxAmplitude;
		
		// writeNodeResults
		std::string nodeSet;
		std::list<std::string> fields;
		std::list<int> timeSteps;
		std::string timeStepTag;
		
		// writeElementResults
		std::string elementSet;
		std::string position;
		
		// writeModalResults
		bool writeModes;
		
		// writeElementProperties
		std::list<std::string> properties;
		
		// writeObjective
		std::list<std::string> objInclude;
		bool writeGradient;
		
		JobCommand();
};

#endif