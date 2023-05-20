#ifndef DVARCLASS
#define DVARCLASS
#include <string>
#include "ListEntClass.h"
#include "DiffDoubClass.h"

class DesignVariable {
	private:
	    std::string category;
		int component;
		int layer;
		double activeTime[2];
		std::string elSetName;
		std::string ndSetName;
		DoubList coefs;
		Doub value;
		DiffDoub diffVal;
		DesignVariable *nextDV;
		
    public:
	    DesignVariable(std::string newCat);
		
		void setActiveTime(double newAt[]);
		
		void setLayer(int newLay);
		
		void setElset(std::string newElset);
		
		void setNdset(std::string newNdset);
		
		void setValue(double newVal);
		
		void setDiffVal(double newVal, double dNewVal);
		
		void setComponent(int newComp);
		
		void setNext(DesignVariable* newNext);
		
		void addCoefficient(double newCoef);		
		
		std::string getCategory();
		
		void getValue(Doub& inp);
		
		void getValue(DiffDoub& inp);
		
		int getComponent();
		
		int getLayer();
		
		DesignVariable* getNext();
		
		void destroy();
};

class DVPt {
	public:
	    DesignVariable *ptr;
		
		DVPt();
};

class DesVarList {
	private:
	    DesignVariable *firstDVar;
		DesignVariable *lastDVar;
		int length;
		
	public:
	    DesVarList();
		
		void addDVar(DesignVariable* newDVar);
		
		int getLength();
		
		DesignVariable* getFirst();
};

#endif