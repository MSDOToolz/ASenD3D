#ifndef JOBCLASS
#define JOBCLASS
#include "ListEntClass.h"
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
		bool nonlinearGeom;
		bool saveSolnHist;
		double simPeriod;
		int solverBandwidth;
		int solverBlockDim;
		std::string solverMethod;
		double staticLoadTime;
		bool thermal;
		double timeStep;
		
		// modal
		std::string type;
		int numModes;
		
		// setSolnToMode
		std::string solnField;
		int mode;
		double maxAmplitude;
		
		// writeNodeResults
		std::string nodeSet;
		StringList fields;
		IntList timeSteps;
		
		// writeElementResults
		std::string elementSet;
		
		// writeModalResults
		bool writeModes;
		
		// writeElementProperties
		StringList properties;
		
		// writeObjective
		StringList objInclude;
		bool writeGradient;
		
		JobCommand *next;
		
		JobCommand();
};

class Job {
	private:
	    JobCommand *firstCmd;
		JobCommand *lastCmd;
		int length;
		
	public:
	    Job();
		
		void addCommand(JobCommand *newCmd);
		
		int getLength();
		
		JobCommand* getFirst();
};

#endif