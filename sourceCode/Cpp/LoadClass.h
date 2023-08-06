#ifndef LOADCLASS
#define LOADCLASS
#include <string>
#include <cmath>
#include "SetClass.h"

class Load {
	private:
	    std::string type;
		double activeTime[2];
		std::string nodeSet;
		Set *ndSetPtr;
		std::string elementSet;
		Set *elSetPtr;
		Set single;
		double load[6];
		double normalDir[3];
		double normTol;
		double center[3];
		double axis[3];
		double angularVel;
		Load *nextLd;
		
	public:
	    Load(std::string newType);
		
		void setActTime(double newAt[]);
		
		void setNodeSet(std::string newSet);
		
		void setNdSetPtr(Set *newSet);
		
		void setElSet(std::string newSet);
		
		void setElSetPtr(Set *newSet);
		
		void setLoad(double newLd[]);
		
		void setNormDir(double newNDir[]);
		
		void setNormTol(double newTol);
		
		void setCenter(double newCent[]);
		
		void setAxis(double newAxis[]);
		
		void setAngVel(double newAngVel);
		
		void setNext(Load* newNext);
		
		// get
		
	    std::string getType();
		
		void getActTime(double actTm[]);
		
		std::string getNodeSet();
		
		Set* getNdSetPtr();
		
		std::string getElSet();
		
		Set* getElSetPtr();
		
		void getLoad(double ldOut[]);
		
		void getNormDir(double ndOut[]);
		
		double getNormTol();
		
		void getCenter(double centOut[]);
		
		void getAxis(double axisOut[]);
		
		double getAngVel();
		
		Load* getNext();
};

class LoadList {
	private:
	    Load *firstLoad;
		Load *lastLoad;
		int length;
		
	public:
	    LoadList();
		
		void addLoad(Load *newLd);
		
		int getLength();
		
		Load* getFirst();

		void destroy();
};

#endif