#ifndef LOADCLASS
#define LOADCLASS
#include <string>
#include <cmath>
#include "SetClass.h"

class Load {
	public:
	    std::string type;
		double activeTime[2];
		std::string nodeSet;
		int ndSetPtr;
		std::string elementSet;
		int elSetPtr;
		double load[6];
		double normalDir[3];
		double normTol;
		double center[3];
		double axis[3];
		double angularVel;
		
	    Load();
		
		void setActTime(double newAt[]);
		
		void setLoad(double newLd[]);
		
		void setNormDir(double newNDir[]);
		
		void setCenter(double newCent[]);
		
		void setAxis(double newAxis[]);
};

#endif